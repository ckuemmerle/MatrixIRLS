function [model, infos] = R3MC_adp(prob,params)%(data_ts, data_ls, model, params)
    %[model, infos] = R3MC(data_ts, data_ls, model, params)
    %
    %
    % Parementers
    % -----------
    % model: Stuctured array with low-rank factors
    %
    % data_ls: data_ls.rows: rows of the given entries
    %        : data_ls.cols: columns of the given entries
    %        : data_ls.entries: values of the given entries
    %        : data_ls.nentries: number of given entries
    %
    % data_ts: 
    %        : 
    %
    % params.tol : Tolerance value used in the stopping criterion
    % params.verbose: Verbosity
    % params.sigma_armijo: The constant factor used in the Armijo linesearch
    % params.ls_maxiter: Maximum of linesearch allowed
    % params.beta_type : 'off' means Gradient descnet,
    %                    'P-R' means Conjugate gradient with Polak-Ribiere+,
    %                    'F-R' means Conjugate gradient with Fletcher-Reeves
    %                    'H-S' means Conjugate gradient with Hestenes-Stiefel+
    % params.linearized_linesearch : true for linearized linesearch,
    %                         : false for adaptive stepsize (arXiv:1209.0430, Equation (13))
    % params.accel_linesearch: true for accelerated linesearch, solving
    %                          degree 2 polynomial approximation
    % params.compute_predictions: To compute test predictions at each iteration
    %                             by randomly choosing a subset of entries
    %
    % Output:
    % -------
    % model: Structured array with final low-rank factors
    % infos: Structured array with additional information
    
    % Refer "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
    % B. Mishra and R. Sepulchre,
    % Technical report, arXiv:1306.2672, 2013.
    % This implementation is due to
    % Bamdev Mishra <b.mishra@ulg.ac.be>, 2013
    
    % Set default parameters
    if ~isfield(params,'beta_type'); params.beta_type = 'P-R'; end % Other choices: 'H-S', 'F-R' and 'off'
    if ~isfield(params,'tol'); params.tol = 1e-5; end % Tolerance on absolute value
    if ~isfield(params,'vtol'); params.vtol = 1e-5; end % Tolerance on relative value
    if ~isfield(params,'verbose'); params.verbose = true; end % Verbosity
    if ~isfield(params,'max_step'); params.max_step = 100; end  % Maximum step size during line-search
    if ~isfield(params,'sigma_armijo'); params.sigma_armijo = 0.5; end % Parameter of the Armijo step
    if ~isfield(params,'ls_maxiter'); params.ls_maxiter = 10; end % Maximum number of line-search steps
    if ~isfield(params,'N0_firstorder'); params.N0_firstorder = 500; end % Maximum number of iterations
    if ~isfield(params,'linearized_linesearch'); params.linearized_linesearch = true; end % Exact linesearch or adapative
    if ~isfield(params,'compute_predictions'); params.compute_predictions = true; end % Compute test error
    if ~isfield(params,'orth_value'); params.orth_value = Inf; end % Measure orthogonality while computing conjugate direction
    if ~isfield(params,'gamma'); params.gamma = 0; end  % Regularization parameter
    if ~isfield(params,'accel_linesearch'); params.accel_linesearch = true; end % To use accelerated lineseach
    
    % User defined options
    tol = params.tol;
    verbose = params.verbose;
    N0 = params.N0_firstorder;
    max_step = params.max_step;
    ls_maxiter = params.ls_maxiter;
    gamma = params.gamma;
    ORTH_VALUE = params.orth_value; % Measure of orthogoanlity used in the CG schemes
    accel_linesearch = params.accel_linesearch;
    linearized_linesearch = params.linearized_linesearch;
    compute_predictions = params.compute_predictions;
    
    if accel_linesearch
        guess_linesearch = @guess_linesearch_urv_accel; % Default linesearch is costly
    else
        guess_linesearch = @guess_linesearch_urv; % This linesearch makes one iteration twice costly
        warning('R3MC:linesearch', 'Going for the default full linearized linesearch. It is costly.\n')
    end
    if isfield(params,'saveiterates') && params.saveiterates == 1
        saveiterates = 1;
    else
        saveiterates = 0;
    end
    
    % Problem dimensions
    model.r = prob.r;
    model.d1 = prob.d1;
    model.d2 = prob.d2;
    r = model.r;
    d1 = model.d1;
    d2 = model.d2;
    data_ls = prob.data_ls;
    
    if ~all(isfield(prob,{'U0_init','R0_init','V0_init'}) == 1)
        warning('R3MC:initialization','Initialization warning: Initialization not provided. Initialized randomly.\n');
        G = randn(d1, r); H = randn(d2, r);
        [Qg Rg] = qr(G, 0);
        [Qh Rh] = qr(H, 0);
        [q1 b q2] = svd(Rg * Rh');
        model.U = Qg * q1;
        model.V = Qh * q2;
        model.R = b;
    else
        model.U = prob.U0_init;
        model.V = prob.V0_init;
        model.R = prob.R0_init;
    end
    
    if ~all(isfield(prob,{'M'}) == 1)
        warning('R3MC:sparse_skeleton','Sparse skeleton warning: Created a sparse skeleton.\n');
        model.M = sparse(data_ls.rows, data_ls.cols, data_ls.entries, d1, d2, data_ls.nentries);
    else
        model.M = prob.M;
    end
    
    
    % Compute cost function
    t_begin0 = tic;
    
    n = data_ls.nentries;
    preds_ls = partXY((model.U * model.R)', model.V', data_ls.rows, data_ls.cols, n)';
    errors = (preds_ls - data_ls.entries);
    
    M = model.M;
    updateSparse(M, errors);
    loss = 0.5*(errors'*errors) + 0.5*gamma*norm(model.R,'fro')^2;%   mean(errors.^2);
    
    tignore0 = tic;
    % If interested in computing recovery
    if compute_predictions
        preds_test = partXY((model.U* model.R)', model.V',data_ts.rows,data_ts.cols, data_ts.nentries)';
        errors_test = preds_test - data_ts.entries;
        cost_test = (errors_test'*errors_test)/data_ts.nentries;
        infos.test_error =  cost_test;
    end
    time_ignore0 = toc(tignore0);
    
    cost = loss;
    
    model.RRt = model.R*model.R';
    model.RtR = model.R'*model.R;
    model.invRRt = eye(r)/model.RRt;
    model.invRtR = eye(r)/ model.RtR;
    
    MV = M * model.V;
    MU = M'* model.U;
    grad.U = MV * model.R';
    grad.V = MU * model.R;
    grad.R = (model.U' * MV) + gamma*model.R;
    grad = egrad2rgrad_urv(model, grad); % Euclidean gradient to Riemannian gradient
    
    
    ip_grad_grad = inner_product_urv(model, grad, grad);
    
    if verbose
        fprintf('[%0.4d] Cost = %5.3e\n', 0, (2/n)*cost);
    end
    
    infos.sqnormgrad = [];
    infos.beta = 0;
    infos.linesearch = [];
    infos.costs = (2/n)*cost;
    j = 0;
    
    dir.U = -grad.U;
    dir.V = -grad.V;
    dir.R = -grad.R;
    alpha = 0;
    
    linesearch_fail = 0;
    ip_grad_dir = inner_product_urv(model, dir, grad);
    
    infos.iter_time = toc(t_begin0) - time_ignore0;
    
    if saveiterates
        X           = cell(1,N0);
    end
    
    for iter = 1:N0
        t_begin = tic; % Begin time.
        
        % Compute first step.
        %   If we don't do an exact line search all the time,
        %   we still need it once for the starting value.
        if linearized_linesearch || iter == 1
            model_new = model;
            stepSize_guess = feval(guess_linesearch, model, errors, dir, data_ls); % An approximate search.
            if stepSize_guess == 0 % We do not trust the linesearch
                stepSize_guess = alpha; % This is the last step;
            end
            alpha = stepSize_guess;
            max_step = alpha;
            sigma_armijo = 1e-4; % Keep a low value or 0 for an exact search
            
        else % Adaptive stepsize
            if j == 0,
                max_step = 2*max_step;
            elseif j > 1,
                max_step = 2*alpha;
            end
            model_new = model;
            alpha = max_step;
            alpha = -sign(ip_grad_dir)*abs(alpha); % Sign correction
            sigma_armijo = params.sigma_armijo; % User specified value
        end
        
        Ui = model.U;
        Vi = model.V;
        Ri = model.R;
        
        model_new.U = uf(Ui + alpha*dir.U);
        model_new.V = uf(Vi + alpha*dir.V);
        model_new.R = Ri + alpha*dir.R;
        
        preds_ls = partXY((model_new.U * model_new.R)', model_new.V', data_ls.rows, data_ls.cols, n)';
        errors = (preds_ls - data_ls.entries);
        cost_new = 0.5*(errors'*errors) + 0.5*gamma*norm(model_new.R,'fro')^2;
        armijo = (cost - cost_new) >= sigma_armijo * abs(alpha * ip_grad_dir);
        
        j = 0;
        while ~armijo,
            j = j + 1;
            alpha = alpha / 2;
            
            model_new.U = uf(Ui + alpha*dir.U);
            model_new.V = uf(Vi + alpha*dir.V);
            model_new.R = Ri + alpha*dir.R;
            
            preds_ls = partXY((model_new.U * model_new.R)', model_new.V', data_ls.rows, data_ls.cols, n)';
            
            errors = (preds_ls - data_ls.entries);
            cost_new = 0.5*(errors'*errors) + 0.5*gamma*norm(model_new.R,'fro')^2;
            armijo = (cost - cost_new) >= sigma_armijo * abs(alpha * ip_grad_dir);
            
            if j > ls_maxiter,
                fprintf('*** Linesearch failed. Quitting...\n');
                linesearch_fail = 1;
                break;
            end
        end
        
        updateSparse(M, errors);
        
        model_new.RRt = model_new.R*model_new.R';
        model_new.RtR = model_new.R'*model_new.R;
        model_new.invRRt = eye(r)/model_new.RRt;
        model_new.invRtR = eye(r)/ model_new.RtR;
        
        
        MV_new = M * model_new.V;
        MU_new = M' * model_new.U;
        
        grad_new.U = MV_new * model_new.R';
        grad_new.V = MU_new * model_new.R;
        grad_new.R =(model_new.U' * MV_new) + gamma*model_new.R;
        grad_new = egrad2rgrad_urv(model_new, grad_new);  % Euclidean gradient to Riemannian gradient
        
        ip_grad_grad_new = inner_product_urv(model_new, grad_new, grad_new);
        
        infos.sqnormgrad = [infos.sqnormgrad; ip_grad_grad_new];
        infos.linesearch = [infos.linesearch; j];
        infos.costs = [infos.costs; (2/n)*cost_new];
        
        tignore = tic;
        % If interested in computing recovery
        if compute_predictions,
            preds_test = partXY((model_new.U*model_new.R)', model_new.V',data_ts.rows,data_ts.cols, data_ts.nentries)';
            errors_test = preds_test - data_ts.entries;
            cost_test = (errors_test'*errors_test)/data_ts.nentries;
            infos.test_error = [infos.test_error; cost_test];
        end
        time_ignore = toc(tignore);
        
        if verbose
            fprintf('[%0.4d] Cost = %7.3e, #extra linesearch = %i, Gradient norm sq. = %7.3e, Step = %7.3e \n', iter, (2/n)*cost_new, j, ip_grad_grad_new, alpha);
        end
        
        if strcmpi(params.beta_type, 'off') % Gradient descent
            beta = 0;
            dir.U = -grad_new.U;
            dir.R = -grad_new.R;
            dir.V = -grad_new.V;
            
        else
            grad_old = vector_transport_urv(model, grad, model_new); % vector transport to the new point
            orth_grads = inner_product_urv(model_new, grad_old, grad_new)/ip_grad_grad_new;
            
            if abs(orth_grads) >= ORTH_VALUE
                beta = 0;
                dir.U = -grad_new.U;
                dir.R = -grad_new.R;
                dir.V = -grad_new.V;
                
            else % Compute the CG modification
                dir = vector_transport_urv(model, dir, model_new); % vector transport to the new point
                
                if strcmp(params.beta_type, 'F-R')  % Fletcher-Reeves
                    beta = ip_grad_grad_new / ip_grad_grad;
                    
                elseif strcmp(params.beta_type, 'P-R')  % Polak-Ribiere+
                    % vector grad(new) - transported grad(current)
                    diff.U = grad_new.U - grad_old.U;
                    diff.R = grad_new.R - grad_old.R;
                    diff.V = grad_new.V - grad_old.V;
                    
                    ip_diff = inner_product_urv(model_new, grad_new, diff);
                    
                    beta = ip_diff / ip_grad_grad;
                    beta = max(0, beta);
                    
                elseif strcmp(params.beta_type, 'H-S')  % Hestenes-Stiefel+
                    diff.U = grad_new.U - grad_old.U;
                    diff.V = grad_new.V - grad_old.V;
                    diff.R = grad_new.R - grad_old.R;
                    
                    ip_diff = inner_product_urv(model_new, grad_new, diff);
                    ip_diff_dir = inner_product_urv(model_new, diff, dir);
                    beta = ip_diff / ip_diff_dir;
                    beta = max(0, beta);
                else
                    error(['Unknown params.beta_type. ' ...
                        'Should be off, F-R, P-R or H-S.']);
                end
                
                dir.U = -grad_new.U + beta*dir.U;
                dir.R = -grad_new.R + beta*dir.R;
                dir.V = -grad_new.V + beta*dir.V;
            end
            
        end
        
        ip_grad_dir_new = inner_product_urv(model_new, dir, grad_new);
        
        % If CG update produces an accent direction, we either reverse it or
        % switch to the negative gradient direction.
        % You came to the wrong neighborhood, bro!
        if sign(ip_grad_dir_new) > 0, %ip_grad_dir_new /ip_grad_grad_new > -1e-10; % sign(ip_grad_dir_new) > 0,
            %         % Either we reverse it.
            %         ip_grad_dir_new = -ip_grad_dir_new;
            %         dir.U = -dir.U;
            %         dir.V = -dir.V;
            %         dir.R = -dir.R;
            %         if verb,
            %             warning('R3MC:accent', 'An accent direction obtained. Reversed it.\n');
            %         end
            
            % Or we switch to the negative gradient direction.
            dir.U = -grad_new.U;
            dir.V = -grad_new.V;
            dir.R = -grad_new.R;
            beta = 0;
            ip_grad_dir_new = -ip_grad_grad_new;
            if verbose
                warning('R3MC:accent','CG warning: An accent direction obtained. Switched to negative gradient dir.\n');
            end
            
        end
        
        % Note: By this time we should have properly chosen the CG direction.
        % Next: We collect the information and prepare for the next iteration.
        
        infos.beta = [infos.beta; beta];
        
        grad = grad_new;
        ip_grad_grad = ip_grad_grad_new;
        ip_grad_dir = ip_grad_dir_new;
        model = model_new;
        infos.iter_time = [infos.iter_time; toc(t_begin) - time_ignore]; % Per iteration time.
        if saveiterates
            X{iter} = {model.U*model.R,model.V};
        end
        if 2*cost_new < tol || linesearch_fail || abs(cost - cost_new)/cost < tol %(2/n)*cost_new < tol
            infos.iter_time = [infos.iter_time; toc(t_begin) - time_ignore]; % Per iteration time.
            N = iter;
            break;
        else
            N = iter;
        end
        cost = cost_new;

    end
    
    model = model_new;
    
    myitertime = infos.iter_time;
    infos.iter_time = cumsum(myitertime);
    if saveiterates
        infos.X = X(1:N);
    else
        infos.X = {model.U*model.R,model.V};
    end
    infos.time = infos.iter_time(1:N);
    infos.N = N;
end


