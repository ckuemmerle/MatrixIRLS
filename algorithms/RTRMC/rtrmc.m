function [U, W, stats] = rtrmc(problem, opts, U0)
%
% Implementation of the Riemannian trust-region methods for matrix
% completion described in [1].
% Code by Nicholas Boumal, cf.
% http://web.math.princeton.edu/~nboumal/RTRMC/index.html
% =========================================================================
% [1] Nicholas Boumal and Pierre-Antoine Absil, Low-rank matrix completion
% via preconditioned optimization on the grassmann manifold. 
% Linear Algebra and its Applications, 15(475): 200?239, 2015.
% =========================================================================
% Modifications by Christian Kuemmerle:
% - save intermediate iterates if opts.saveiterates == 1. An according
% modification has also been done 
% =========================================================================
% Given a problem structure describing a low-rank matrix completion problem
% instance, as obtained from the buildproblem function, will return U and
% W, such that U is orthonormal and the product U*W is an estimate of the
% sought matrix. The rank of the approximation U*W is given by problem.r.
% The opts structure is optional: see in code for available options.
% U0 is an initial guess of an orthonormal matrix whose columns span the
% column space of the sought matrix. If not specified or left empty, this
% will be created using the initialguess function.
%
% Algorithm timings should be read from the stats structure returned, as it
% makes the distinction between time actually spent in the algorithm and
% time spent computing debugging quantities such as the RMSE at each
% iteration.
%
% By default, this function behaves like RTRMC 2p. Play with the options
% opts.method, opts.order and opts.precon to obtain the other algorithms:
%  RTRMC 2  : method = 'rtr', order = 2, precon = false
%  RTRMC 2p : method = 'rtr', order = 2, precon = true
%  RTRMC 1  : method = 'rtr', order = 1, precon = false
%  RCGMC    : method = 'cg',  order = 1, precon = false
%  RCGMCp   : method = 'cg',  order = 1, precon = true
%
% Algorithm by Nicolas Boumal and P.-A. Absil
% Code by Nicolas Boumal, UCLouvain, 2013.
%
% See also: buildproblem initialguess

    if exist('manoptsolve', 'file') ~= 2
        error(['Please add the Manopt toolbox to your path first. ' ...
               'You can do so by calling the importmanopt script in ' ...
               'the Manopt_2.0 folder packaged with this release, ' ...
               'or by installing the toolbox from www.manopt.org.']);
    end
    
    % Gather values from the problem description
    m = problem.m;
    n = problem.n;
    r = problem.r;
    k = problem.k;
	X = problem.X;
    C = problem.C;
    I = problem.I;
    J = problem.J;
    mask = problem.mask;
	Chat = problem.Chat;
	lambda = problem.lambda;

    % Pick the standard initial guess if none is given.
    if nargin < 3 || isempty(U0)
        tic_initguess = tic();
        U0 = initialguess(problem);
        time_initguess = toc(tic_initguess);
    else
        time_initguess = 0;
    end
    
    % Define the default options of the method.
    default_opts = struct('method', 'rtr', ...
                          'order', 2, ...
                          'precon', true, ...
                          ... %
                          ... % See Manopt's documentation for more options
                          ... %
                          'maxiter', 300, ...
                          'maxinner', 500, ...
                          'mininner', 0, ...
                          'tolgradnorm', 1e-5, ...
                          'useRand', false, ...
                          'kappa', 0.1, ...
                          'theta', 1, ...
                          'rho_prime', 0.1, ...
                          'storedepth', 5, ...
                          'computeRMSE', false, ...
                          'verbosity', 2);
    
    % Merge them with the user supplied option values, if any. The default
    % options will be complemented with (and possibly overwrittin by) the
    % user supplied options.
    if ~exist('opts', 'var') || isempty(opts)
        opts = struct();
    end
    opts = mergeOptions(default_opts, opts);
    
    % For compatibility with legacy test code ...
    if isfield(opts, 'maxouter')
        opts.maxiter = opts.maxouter;
    end
    
    % For RTRMC 1, we disable rho regularization, as it prevents fine
    % convergence. This deserves further investigation.
    if ~isfield(opts, 'rho_regularization')
        if strcmpi(opts.method, 'rtr') && opts.order == 1
            opts.rho_regularization = 0;
        end
    end
    
    % Setup the initial trust region radius and the maximum trust region
    % radius depending on the preconditioner. The preconditioner changes
    % the shape of the trust region, and hence also changes the lengths of
    % the vectors which are to be compared to Delta.
    % Mnorm is the 2-norm (the largest eigenvalue) of the preconditioner
    % before it is inverted.
    if opts.precon
        % We do not account for the time it takes to computes W0 here,
        % since it will be computed in the optimization algorithm anyway,
        % and timed there.
        W0 = lsqfit(problem, U0);
        Mnorm = sqrt(norm((W0*W0'), 2)/k);
    else
        Mnorm = 1;
    end
    if ~isfield(opts, 'Delta_bar')
        % Base length: the diameter of the Grassmann manifold.
        opts.Delta_bar = Mnorm * sqrt(r)*pi/2;
    end
    if ~isfield(opts, 'Delta0')
        opts.Delta0 = opts.Delta_bar / 8;
    end

    % Check some of the option values
    if opts.order ~= 1 && opts.order ~= 2
        warning('RTRMC:order', ...
          ['The RTRMC ''order'' option must be set to 1 (no Hessian) ' ...
           'or 2 (with Hessian).']);
        opts.order = 2;
    end
    
    if m > n
        warning(['RTRMC is optimized for m <= n (matrices with more ' ...
                 'columns than rows). Please transpose the problem.']); %#ok<WNTAG>
    end
    
    % Obtain a description of the Grassmann manifold for Manopt
    grass = grassmannfactory(m, r);
    optproblem.M = grass;
    
    % Specify the cost and the gradient functions (see below)
    optproblem.cost = @cost;
    optproblem.grad = @gradient;

    % Specify whether to use the Hessian or not. If not, "approximate" it
    % with the identity matrix.
    if opts.order == 2
        optproblem.hess = @hessian;
    elseif opts.order == 1
        optproblem.approxhess = @(U, H) H;
    end
    
    % Preconditioner for the Hessian.
    % Be careful: if W becomes rank-deficient or close, this will crash.
    % A rank-deficient W is the sign that r should be lowered. This is not
    % currently monitored in the algorithm, but could easily be.
    if opts.precon
        optproblem.precon = @precon;
        optproblem.sqrtprecon = @sqrtprecon;
    end

    % RTRMC can compute the RMSE at each iteration efficiently if required,
    % either for the full matrix X=AB if the factors A and B are provided,
    % or some test RMSE if data for that is provided (see recordrmse).
    % Time spent in statsfun is not included in the algorithm timing.
    if opts.computeRMSE
        opts.statsfun = @recordrmse;
    end
    
    if opts.saveiterates
       opts.recorditerates = @recorditerates;
    end
    
  
    % Do the magic.
    if strcmpi('rtr', opts.method)
        [U, bestcost, stats] = trustregions(optproblem, U0, opts); %#ok<ASGLU>
    elseif strcmpi('cg', opts.method)
        if ~isfield(opts, 'beta_type')
            opts.beta_type = 'H-S';
        end
        opts.linesearch = @linesearch;
        opts.ls_contraction_factor = .2;
        opts.ls_optimism = 1.1;
        opts.ls_suff_decr = 1e-4;
        opts.ls_max_steps = 25;
        [U, bestcost, stats] = conjugategradient(optproblem, U0, opts); %#ok<ASGLU>
    end

    % Account for the initial guess computation time
    for i = 1 : length(stats)
        stats(i).time = stats(i).time + time_initguess;
    end
    
    % Return also the W factor of the factorization X ~ UW.
    % We do not account for the time it takes to computes W here, since it
    % was computed in the optimization algorithm anyway, and timed there.
    W = lsqfit(problem, U);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Define the cost function
    function [f, store] = cost(U, store)
        
        if ~isfield(store, 'val')
            
            [W, store] = lsqfit(problem, U, store);
            store.W = W;
            
            UW = spmaskmult(U, W, I, J);
            store.UW = UW;
            
            store.val = ( .5*sum((C.*(UW-X)).^2)    ...
                        + .5*lambda.^2*sum(W(:).^2) ...
                        - .5*lambda.^2*sum(UW.^2)   ) / k;
        end
        f = store.val;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define the gradient function
    function [grad, store] = gradient(U, store)
        
        if ~isfield(store, 'grad')

            % If the cost was never computed for this point before, call it
            % so we have the necessary ingredients to compute the gradient.
            if ~isfield(store, 'val')
                [~, store] = cost(U, store);
            end
        
            W = store.W;
            WWt = W*W.';
            sqlaWWt = lambda.^2*WWt;
            store.WWt = WWt;
            store.sqlaWWt = sqlaWWt;

            UW = store.UW;
            RU = Chat.*(UW-X) - (lambda.^2)*X;
            store.RU = RU;
            
            store.grad = (multsparsefull(RU, W.', mask) + U*sqlaWWt)/k;
            
        end
        grad = store.grad;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define the Hessian function
    function [hess, store] = hessian(U, H, store)
        
        [W, dW, store] = lsqfit(problem, U, H, store);
        UdW = spmaskmult(U, dW, I, J);
        
        if ~isfield(store, 'HW')
            store.HW  = spmaskmult(H, W, I, J);
        end
        HW = store.HW;
        
        if ~isfield(store, 'grad')
            [~, store] = gradient(U, store);
        end
        RU = store.RU;
        sqlaWWt = store.sqlaWWt;

        hess = multsparsefull(Chat.*(HW + UdW), W.', mask);
        hess = hess - U*(U.'*hess);
        hess = hess + multsparsefull(RU, dW.', mask);
        hess = hess + H*sqlaWWt;
        hess = hess + U*(lambda.^2*(W*dW.'));
        hess = hess/k;
            
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define the preconditioner
    function [pH, store] = precon(U, H, store)
        
        if ~isfield(store, 'WWt')
            [W, store] = lsqfit(problem, U, store);
            store.W = W;
            WWt = W*W.';
            store.WWt = WWt;
        end
        WWt = store.WWt;
        
        pH = H / (WWt/k);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% The square root of the preconditioner is only useful to compute the
    %  spectrum of the Hessian with and without preconditioning. See
    %  Manopt's hessianspectrum and hessianextreme functions if you want
	%  to compute that.
    function [sqrtpH, store] = sqrtprecon(U, H, store)
        
        if ~isfield(store, 'WWt')
            [W, store] = lsqfit(problem, U, store);
            store.W = W;
            WWt = W*W.';
            store.WWt = WWt;
        end
        WWt = store.WWt;
        
        sqrtpH = H / (sqrtm(WWt)/sqrt(k));
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Keep track of the RMSE
    function stats = recordrmse(optproblem, U, stats, store) %#ok<INUSL>
        if ~isfield(store, 'W')
            W = lsqfit(problem, U);
        else
            W = store.W;
        end
        [rmse, rmse_clipped] = computeRMSE(U, W, problem);
        stats.RMSE = rmse;
        stats.RMSE_clipped = rmse_clipped;
        stats.WWt_cond = cond(store.WWt);
        if opts.verbosity > 2
            fprintf('cond(W*W'') : %e\n', stats.WWt_cond);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Keep track of the RMSE
    function stats = recorditerates(optproblem, U, stats, store) %#ok<INUSL>
        stats.U=U;
        stats.W=store.W;
%         if ~isfield(store, 'W')
%             W = lsqfit(problem, U);
%         else
%             W = store.W;
%         end
%         [rmse, rmse_clipped] = computeRMSE(U, W, problem);
%         stats.RMSE = rmse;
%         stats.RMSE_clipped = rmse_clipped;
%         stats.WWt_cond = cond(store.WWt);
%         if opts.verbosity > 2
%             fprintf('cond(W*W'') : %e\n', stats.WWt_cond);
%         end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [rmse, rmse_clipped] = computeRMSE(U, W, problem)

        if isfield(problem, 'A') && isfield(problem, 'B')

            A = problem.A;
            B = problem.B;
            m = problem.m;
            n = problem.n;
            rmse = sqrt(sqfrobnormfactors(U, W, A, B)/(m*n));
            rmse_clipped = rmse;

        elseif isfield(problem, 'Xtest') && isfield(problem, 'Itest') && ...
               isfield(problem, 'Xmean') && isfield(problem, 'Jtest')

            Xtest = problem.Xtest;
            Itest = problem.Itest;
            Jtest = problem.Jtest;
            Xmean = problem.Xmean;
            Xpred = Xmean + spmaskmult(U, W, Itest, Jtest);
            Xpred_clipped = max(1, min(5, Xpred));
            rmse = norm(Xpred - Xtest, 'fro')/sqrt(length(Xtest));
            rmse_clipped = norm(Xpred_clipped - Xtest, 'fro')/sqrt(length(Xtest));

        else
            rmse = NaN;
            rmse_clipped = NaN;
        end

        if opts.verbosity > 2
            fprintf('Test RMSE : %e\n', rmse);
            fprintf('Test RMSE clipped : %e\n', rmse_clipped);
        end

    end

end
