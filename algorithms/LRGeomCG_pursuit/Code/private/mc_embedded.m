function [U stats optproblem] = mc_embedded(problem, opts, X0, solver)
% Solver field:
%   'cg' = non-linear CG
%   'rtr-exact-hessian' = RTR with exact Riemannian Hessian
%   'rtr-fd-hessian'    = RTR with finite difference Hessian
%   'rtr-ehess2rhess'   = RTR with exact Hessian obtained via ehess2rhess

    if exist('manopt_version', 'file') ~= 2
        error(['Please add the Manopt toolbox to your path first. ' ...
               'You can do so by calling the importmanopt script in ' ...
               'the Manopt_1.0.5 folder packaged with this release, ' ...
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
    time_initguess = 0;
    
    
    % Define the default options of the method.
    default_opts = struct('method', 'cg', ...
                          'order', 2, ...
                          'precon', false, ...
                          ... %
                          ... % See Manopt's documentation for more options
                          ... %
                          'maxiter', 1000, ...
                          'maxinner', 500, ...
                          'mininner', 0, ...
                          'tolgradnorm', 1e-5, ...
                          'useRand', false, ...
                          'kappa', 0.1, ...
                          'theta', 1, ...
                          'rho_prime', 0.1, ...
                          'storedepth', 5, ...
                          'computeRMSE', true, ...
                          'verbosity', 1);
    
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
    
    % Setup the initial trust region radius and the maximum trust region
    % radius depending on the preconditioner. The preconditioner changes
    % the shape of the trust region, and hence also changes the lengths of
    % the vectors which are to be compared to Delta.
    % Mnorm is the 2-norm (the largest eigenvalue) of the preconditioner
    % before it is inverted.
    if opts.precon
        error('not implemented')
        % We do not account for the time it takes to computes W0 here,
        % since it will be computed in the optimization algorithm anyway,
        % and timed there.
        W0 = lsqfit(problem, X0);
        Mnorm = sqrt(norm((W0*W0'), 2)/k);
    else
        Mnorm = 1;
    end
%     if ~isfield(opts, 'Delta_bar')
%         % Base length: the diameter of the Grassmann manifold.
%         opts.Delta_bar = Mnorm * sqrt(r)*pi/2;
%     end
%     if ~isfield(opts, 'Delta0')
%         opts.Delta0 = opts.Delta_bar / 8;
%     end

    % Check some of the option values
    if opts.order ~= 1 && opts.order ~= 2
        warning('RTRMC:order', ...
          ['The RTRMC ''order'' option must be set to 1 (no Hessian) ' ...
           'or 2 (with Hessian).']);
        opts.order = 2;
    end
    
    
    
    % Obtain a description of the fixedrank manifold for Manopt
    fr_manifold = fixedrankembeddedfactory(m, n, r);    
    optproblem.M = fr_manifold;
    
    % Specify the cost and the gradient functions (see below)
    optproblem.cost = @cost;
    optproblem.grad = @Riemannian_gradient;
    optproblem.egrad = @Euclidean_gradient;
    
%%% FOR INCLUSION IN MANTOPT
    optproblem.M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, H)
        % egrad is a possibly sparse matrix of the same matrix size as X
        % ehess is the matvec of Euclidean Hessian with H, is also a matrix
        % of the same matrix size as X

        rhess = fr_manifold.proj(X, ehess);   % probably other name here...                             

        % Curvature part            
        T = (egrad*H.Vp)/X.S;
        rhess.Up = rhess.Up + (T - X.U*(X.U'*T));

        T = (egrad'*H.Up)/X.S;
        rhess.Vp = rhess.Vp + (T - X.V*(X.V'*T));
    end

    optproblem.M.tangent2ambient = @tangent2ambient;
    function Z = tangent2ambient(x, u)
        Z = x.U*u.M*x.V' + u.Up*x.V' + x.U*u.Vp';
    end


    
    if strcmpi(solver, 'rtr-exact-hessian')
        opts.method = 'rtr';
        optproblem.hess = @Riemannian_hessian;
    elseif strcmpi(solver, 'rtr-fd-hessian')
        opts.method = 'rtr';
        % deliberately no Hessian supplied        
    elseif strcmpi(solver, 'ehess2rhess')    
        opts.method = 'rtr';
        optproblem.ehess = @Euclidean_hessian;
        
        % todo
    elseif strcmpi(solver, 'rtr-only-Euclidean-hessian')
        opts.method = 'rtr';
        optproblem.hess = @proj_Euclidean_hessian; % deliberately "wrong" Hessian supplied        
    elseif strcmpi(solver, 'cg') 
        opts.method = 'cg';        
    elseif strcmpi(solver, 'newton-exact-hessian')
        opts.method = 'newton';
        %optproblem.hess = @Riemannian_hessian;
        optproblem.hess = @proj_Euclidean_hessian; % deliberately "wrong" Hessian supplied        
    else
        error('Set solver option correctly')
    end
    
    
    
    % Preconditioner for the Hessian.
    % Be careful: if W becomes rank-deficient or close, this will crash :/.
    % A rank-deficient W is the sign that r should be lowered. This is not
    % currently monitored in the algorithm, but could easily be.
%     if opts.precon
%         optproblem.precon = @precon;
%     end

    % RTRMC can compute the RMSE at each iteration efficiently if required,
    % either for the full matrix X=AB if the factors A and B are provided,
    % or some test RMSE if data for that is provided (see recordrmse).
    % Time spent in statsfun is not included in the algorithm timing.
    if opts.computeRMSE
        opts.statsfun = @recordrmse;
    end
    
    
    
    
    % Do the magic.
    if strcmpi('rtr', opts.method)
        [U, bestcost, stats] = trustregions(optproblem, X0, opts); %#ok<ASGLU>
    elseif strcmpi('cg', opts.method)
        if ~isfield(opts, 'beta_type')
            opts.beta_type = 'P-R';            
        end
        opts.linesearch = @linearized_linesearch;
        opts.ls_contraction_factor = .5;
        opts.ls_optimism = 1.1;
        opts.ls_suff_decr = 0.001;
        opts.ls_max_steps = 25;
        [U, bestcost, stats] = conjugategradient(optproblem, X0, opts); %#ok<ASGLU>
        %[U, bestcost, stats] = steepestdescent(optproblem, U0, opts); %#ok<ASGLU>
    elseif strcmpi('newton', opts.method)
        opts.linesearch = @linearized_linesearch;
        opts.ls_contraction_factor = .5;
        opts.ls_optimism = 1.1;
        opts.ls_suff_decr = 0.001;
        opts.ls_max_steps = 25;
        [U, bestcost, stats] = newton(optproblem, X0, opts);
    end

    % Account for the initial guess computation time
    for i = 1 : length(stats)
        stats(i).time = stats(i).time + time_initguess;
    end
    
%     t=tic();
%    lam = hessianspectrum(optproblem, U); %, sqrtprec)
%     toc(t)
%     
%     t=tic();
%    lam2 = hessianspectrum_full(optproblem, U); %, sqrtprec)
%     toc(t)
%     
    
    
%  figure(3)
%  hist(log10(lam))
%  pause
% %checkgradient(optproblem, U0);
% checkhessian(optproblem, U);
% pause

% if strcmpi(solver, 'rtr-exact-hessian')
% if strcmpi(solver, 'rtr-only-Euclidean-hessian')
%     figure(10)
%     checkgradient(optproblem, X0);
%     figure(11)
%     checkhessian(optproblem, X0);
%     pause
% end


   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Define the cost function
    function [f store] = cost(x, store)
        
        if ~isfield(store, 'val')
                        
            UW = spmaskmult(x.U*x.S, x.V', I, J); % sampled matrix U
            store.UW = UW;
            RU = (UW-X);
            store.RU = RU(:);
            
            store.val = 0.5*(RU'*RU);
            
            % ugly hack for the linearized_linesearch: we need I and J
            store.I = I; store.J = J;
        end
        f = store.val;
                
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define the Riemannian gradient function
    function [grad store] = Riemannian_gradient(x, store)
        
        if ~isfield(store, 'grad')

            % If the cost was never computed for this point before, call it
            % so we have the necessary ingredients to compute the gradient.
            if ~isfield(store, 'val')
                [~, store] = cost(x, store);
            end                    

            
            setsparseentries(mask, store.RU);                         
            
            store.grad = fr_manifold.egrad2rgrad(x, mask);
            
        end
        grad = store.grad;
        
    end

    %% Define the Riemannian gradient function
    function [egrad store] = Euclidean_gradient(x, store)
        
        if ~isfield(store, 'egrad')

            % If the cost was never computed for this point before, call it
            % so we have the necessary ingredients to compute the gradient.
            if ~isfield(store, 'val')
                [~, store] = cost(x, store);
            end                    

            
            setsparseentries(mask, store.RU);                         
            
            store.egrad = mask*1.0; % copy mask
            
        end
        egrad = store.egrad;
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define the Hessian function
    
    % Riemannian Hessian directly
    function [hess store] = Riemannian_hessian(x, d, store)
        if ~isfield(store, 'val')
            UW = spmaskmult(x.U*x.S, x.V', I, J); % sampled matrix U
            store.UW = UW;
            RU = (UW-X);
            store.RU = RU(:);
            store.val = 0.5*(RU'*RU);
        end
                

        % Euclidean part
        D_omega = spmaskmult([x.U*d.M+d.Up x.U],[x.V d.Vp]', I, J); 
        setsparseentries(mask, D_omega);          
        hess = fr_manifold.proj(x, mask);
                
        setsparseentries(mask, store.RU);  
        
        % Curvature part            
        T = (mask*d.Vp)/x.S;
        hess.Up = hess.Up + (T - x.U*(x.U'*T));
                                
        T = (mask'*d.Up)/x.S;
        hess.Vp = hess.Vp + (T - x.V*(x.V'*T));
    end

    % only Euclidean Hessian directly
    function [ehess store] = Euclidean_hessian(x, d, store)        
        % Euclidean part
        D_omega = spmaskmult([x.U*d.M+d.Up x.U],[x.V d.Vp]', I, J); 
        setsparseentries(mask, D_omega);          
        ehess = mask*1.0; % we need to copy the mask
    end

    % only Euclidean Hessian directly
    function [ehess store] = proj_Euclidean_hessian(x, d, store)        
        % Euclidean part
        D_omega = spmaskmult([x.U*d.M+d.Up x.U],[x.V d.Vp]', I, J); 
        setsparseentries(mask, D_omega);          
        ehess = fr_manifold.proj(x, mask);   
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define the preconditioner
    function [pH store] = precon(U, H, store)
        error('not implemented')
        
        if ~isfield(store, 'WWt')
            [W store] = lsqfit(problem, U, store);
            store.W = W;
            WWt = W*W.';
            store.WWt = WWt;
        end
        WWt = store.WWt;
        
        pH = H / (WWt/k);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% The square root of the preconditioner is only useful to compute the
    %% spectrum of the Hessian with and without preconditioning. See
    %% Manopt's hessianspectrum function if you want to compute that.
    function [sqrtpH store] = sqrt_precon(U, H, store) %#ok<DEFNU>
        error('not implemented')
        
        if ~isfield(store, 'WWt')
            [W store] = lsqfit(problem, U, store);
            store.W = W;
            WWt = W*W.';
            store.WWt = WWt;
        end
        WWt = store.WWt;
        
        sqrtpH = H / (sqrtm(WWt)/sqrt(k));
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Keep track of the RMSE
    function stats = recordrmse(optproblem, x, stats, store) %#ok<INUSL>
        
        [rmse, rmse_clipped] = computeRMSE(x, problem);
        stats.RMSE = rmse;
        stats.RMSE_clipped = rmse_clipped;                
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [rmse rmse_clipped] = computeRMSE(x, problem)

        if isfield(problem, 'A') && isfield(problem, 'B')
            error('not implemented')
            A = problem.A;
            B = problem.B;
            m = problem.m;
            n = problem.n;
            rmse = sqrt(sqfrobnormfactors(x, W, A, B)/(m*n));
            rmse_clipped = rmse;

        elseif isfield(problem, 'Xtest') && isfield(problem, 'Itest') && ...
               isfield(problem, 'Xmean') && isfield(problem, 'Jtest')

            Xtest = problem.Xtest;
            Itest = problem.Itest;
            Jtest = problem.Jtest;
            Xmean = problem.Xmean;
            Xpred = Xmean + spmaskmult(x.U*x.S, x.V', Itest, Jtest);
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




    function [stepsize newx storedb lsmem lsstats] = linearized_linesearch(problem, x, d, f0, df0, opts, storedb, lsmem)                       
        % Adapted from manopt's standard linesearch.m with added functionality:
        %   1) take initial step as the exact minimizer of the linearized version 
        %      (that is, only on the tangent space).
        %   2) force the first step as the suggested step from above.

        % new code (next 5 lines)
        store = getStore(problem, x, storedb);       
        dir_omega = spmaskmult([x.U*d.M+d.Up x.U],[x.V d.Vp]', I, J); 
        A1 = dir_omega;
        tinit = -A1 \ store.RU(:);
        options.ls_initial_stepsize = tinit;




        % Backtracking default parameters. These can be overwritten in the
        % options structure which is passed to the solver.
        default_options.ls_contraction_factor = .5;
        default_options.ls_optimism = 1/.5;
        default_options.ls_suff_decr = 1e-4;
        default_options.ls_max_steps = 25;
        default_options.ls_initial_stepsize = 1;

        options = mergeOptions(default_options, options);

        contraction_factor = options.ls_contraction_factor;
        optimism = options.ls_optimism;
        suff_decr = options.ls_suff_decr;
        max_ls_steps = options.ls_max_steps;
        initial_stepsize = options.ls_initial_stepsize;

        % Compute the norm of the search direction.
        % This is useful to make the linesearch algorithm invariant under the
        % scaling of d. The rationale is that the important information is the
        % search direction, not the size of that vector. The question of how
        % far we should go is precisely what the linesearch algorithm is
        % supposed to answer: the calling algorithm should not need to care.
        norm_d = problem.M.norm(x, d);

        % new code: just use the initial guess and disable the fancy original stuff
        alpha = tinit;

        %     % If we know about what happened at the previous step, we can leverage
        %     % that to compute an initial guess for the step size, as inspired from
        %     % Nocedal&Wright, p59.
        %     if isstruct(lsmem) && isfield(lsmem, 'f0')
        %         % Pick initial step size based on where we were last time,
        %         alpha = 2*(f0-lsmem.f0)/df0;
        %         % and go look a little further (or less far), just in case.
        %         alpha = optimism*alpha;
        %     end
        %     
        %     % If we have no information about the previous iteration (maybe this is
        %     % the first one?) or if the above formula gave a too small step size
        %     % (perhaps it is even negative), then fall back to a user supplied
        %     % suggestion for the first step size (the "a priori").
        %     % At any rate, the choice should be invariant under rescaling of the
        %     % cost function f and of the search direction d, and it should be
        %     % bounded away from zero for convergence guarantees. We must allow it
        %     % to be close to zero though, for fine convergence.
        %     if isnan(alpha) || alpha*norm_d <= eps
        %         alpha = initial_stepsize/norm_d;
        %     end
        %     
        %     if isfield(options, 'force_initial_step') && options.force_initial_step 
        %         alpha = initial_stepsize;
        %     end


        % Make the chosen step and compute the cost there.
        newx = problem.M.retr(x, d, alpha);
        [newf storedb] = getCost(problem, newx, storedb);
        cost_evaluations = 1;

        % Backtrack while the Armijo criterion is not satisfied
        while newf > f0 + suff_decr*alpha*df0

            % Reduce the step size,
            alpha = contraction_factor * alpha;

            % and look closer down the line
            newx = problem.M.retr(x, d, alpha);
            [newf storedb] = getCost(problem, newx, storedb);
            cost_evaluations = cost_evaluations + 1;

            % Make sure we don't run out of budget
            if cost_evaluations >= max_ls_steps
                break;
            end

        end

        % If we got here without obtaining a decrease, we reject the step.
        if newf > f0
            alpha = 0;
            newx = x;
            newf = f0; %#ok<NASGU>
        end

        % As seen outside this function, stepsize is the size of the vector we
        % retract to make the step from x to newx. Since the step is alpha*d:
        stepsize = alpha * norm_d;

        % Save the situtation faced now so that at the next iteration, we will
        % know something about the previous decision.
        lsmem.f0 = f0;
        lsmem.df0 = df0;
        lsmem.stepsize = stepsize;

        % Save some statistics also, for possible analysis.
        lsstats.costevals = cost_evaluations;
        lsstats.stepsize = stepsize;
        lsstats.alpha = alpha;
    
    end

end
