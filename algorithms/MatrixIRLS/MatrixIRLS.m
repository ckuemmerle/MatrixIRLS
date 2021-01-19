function [Xr,outs] = ...
    MatrixIRLS(prob,lambda,opts)
% Given partial observations of a matrix, attempts to find a low-rank completion.
%
% Implements the algorithm Matrix Iteratively Reweighted Least Squares
% (MatrixIRLS) for optimizing a smoothed log-det/Schatten-p objective with
% updated smoothing, as described in the paper
% [1] C. Kuemmerle, C. M. Verdun, "Escaping Saddle Points in 
% Ill-Conditioned Matrix Completion with a Scalable Second Order Method", 
% ICML 2020 Workshop "Beyond First Order Methods in ML Systems".
%
% See also: Chapter 2 of 
% [2] C. Kuemmerle, "Understanding and Enhancing Data Recovery Algorithms From
% Noise-Blind Sparse Recovery to Reweighted Methods for Low-Rank Matrix 
% Optimization", Ph.D. dissertation, Technical University of Munich,
% https://mediatum.ub.tum.de/doc/1521436/1521436.pdf.
%
% In fact, this implementation corresponds to the very algorithm of the
% paper [1] if the second input parameter fulfills lambda == 0. For lambda
% > 0, it implements a variant of MatrixIRLS that does not involves an exact
% data constraint, but rather the optimization of the unconstrained objective
%       F_{\epsilon}(X) + \frac{1}{2 \lambda} \|P_{\Omega}(X)-y\|_2^2.  (1)
% =========================================================================
% Parameters
% ----------
% prob: struct. Contains problem parameters, fields include:
%       - d1        int. First dimension of matrix to completed
%       - d2        int. Second dimension of matrix to completed
%       - Phi       (d1 x d2) sparse matrix. Non-zero entries
%                       correspond to observed matrix entries.
%       - y         (m x 1) vector. Provided matrix entry values,
%                       ordered in according to linear indices of non-zero
%                       entries of prob.Phi.
%       - r         int. (Target) rank/ rank estimate \tilde{r}, to be used 
%                       in the definition of smoothing parameter update 
%                       (if opts.mode_eps=='oracle_model_order'). 
%                       Can be also left empty.
% lambda: double. Data fit parameter, choose as 0 for equality data
%           constraint, and larger value for noisy/inaccurate observations. 
%           Large lambda corresponds to a loose data fit, see (1) above.
% opts: struct. Contains algorithmic options. Fields include:
%       - p         double, between 0 and 1. Parameter of (non-)convexity of
%                   objective function.
%                   = 0: sum of log objective /log-det objective.
%                   > 0 and < 1: Schatten-p quasi-norm objective.
%                   = 1: Nuclear norm type objective.
%       - N0        int. Max. number of (outer) IRLS iterations 
%                   (default 200). 
%       - N0_inner  int. Max. number of inner iteration in the conjugate
%                   gradient method (default 200).
%       - tol       double. Stopping criterion, lower bound for relative 
%                   change of Frobenius norm (default 1e-9).
%       - tol_CG_fac  double. Factor that determines the stopping criterion 
%                   if the conjugate gradient method, such that
%                   tol_CG_fac*eps_c is a lower bound on the relative
%                   residual (eps_c is the current smoothing parameter 
%                   epsilon). (default 1e-5).
%       - epsmin    Minimal value for regularization parameter epsilon
%                   (default 1e-16).
%       - mode_eps   character string. Choice of rule to update the
%                    regularization parameter epsilon (default
%                    'oracle_model_order').
%                   = 'oracle_model_order': Uses knowledge about the
%                      model order/rank "r", epsilon choice as in [1],[2].
%                   = 'auto_decay': eps that is automatically
%                      superlinearly decaying. First value eps0 of epsilon 
%                      is based on the Frobenius norm of first iterate
%                      X^{(l)}, choose eps at iteration l such that:
%                      eps(l) = eps0.*(tol_outer/eps0)^(((l-1)/(N0-1))^(2-p));
%                   = 'iter_diff': epsilon choice rule similar [3], 
%                      based on the norm of the difference of two conse-
%                      cutive iterates.
%       - use_min   Boolean (default 1).
%                   = 1: Force monotonicity of regularization parameter 
%                   epsilon.
%                   = 0: No forced monotonicity, plain application of
%                   epsilon choice rule.
%       - R         int. Maximum rank parameter for intermediate iterates.
%                   (default: min(d1,d2)). Can be set to smaller values 
%                   (upper bound for the rank of the solution), used to 
%                   accelerate algorithm for large problems if target rank 
%                   is misspecified or if a non-oracle choice of the 
%                   smoothing parameter epsilon is chosen.
%       - type_mean   character string, determines type of weight operator
%                   underlying the IRLS method, see [2].
%                   = 'optimal': Corresponds to the optimal weight
%                      operator: Geometric mean of p = 0, and
%                      (p/(p-2))-power mean for Schatten-p objective.
%                   = 'geometric': Geometric mean.
%                   = 'harmonic': Harmonic mean.
%                   = 'arithmetic: Arithmetic mean.
%                   = 'min': Based on minimum of right- and
%                      left-weights.
%                   = 'qmean': q-power mean.
%       - qmean_para  double. Only relevant if type_mean=='qmean', 
%                      determines the value of "q" in the q-power mean,
%                      cf. [2].
%       - increase_antisymmetricweights
%                   Boolean (default = 0). Can be set to 1, but that is
%                   only interesting for theoretical reasons.
%                   = 1: In the definition of the weight operator, this 
%                      means that the weights on the "antisymmetric part" 
%                      are increased as in [2]. 
%                   = 0: Weight operator without artificially increased
%                      antisymmettic port.
%       - objective character string, decides precise form of objective
%                       function to be optimized by IRLS.
%                   = 'objective_thesis': log(max(sigma,eps))+0.5.*
%                       (min(sigma,eps)^2/eps^2-1).
%                       (as used in [1], [2])
%                   = 'pluseps': Corresponds to objective
%                       of the type:
%                       log(max(sigma,eps)+eps)+0.25.*(min(sigma,eps)^2/eps^2-1)
%                   = 'pluseps_squared_max': Corresponds to objectives
%                       of the type:
%                       log(max(sigma,eps)^2+eps^2)+0.5.*(min(sigma,eps)^2/eps^2-1)
%                   = 'pluseps_squared': Corresponds to objectives
%                       of the type:
%                       log(sigma^2+eps^2)  
%       - verbose   int. Detemines the level of displayed text and
%                   amount of output variables (default 1)
%                   = 0: Very little text output, only output basic
%                        parameters such as regularization parameters
%                        eps_c.
%                   = 1: Return also intermediate iterates and a few 
%                        algorithmic parameters, more text display.
%                   = 2: Return also many further parameters.
%       - saveiterates  Boolean. Save or do not save intermediate iterates
%                   (default 1).
%                   = 0: Only output last iterate.
%                   = 1: Output all intermediate iterates.
%       - mode_linsolve  character string, describes linear system to be solved.
%                   = 'tangspace': Linear system based on tangent space,
%                   see Appendix of [1] or Section 2.6.2. of [2].
% Returns
% ----------
% Xr:   struct. Contains the last iterate in space-efficient form, see
%       cofficient matrices in Appendix of [1]. In particular,
%                   Xr.Gam1 = \tilde{Gamma}_1  [(r_k x r_k) matrix]
%                   Xr.Gam2 = \tilde{Gamma}_2  [(r_k x d_2) matrix]
%                   Xr.Gam3 = \tilde{Gamma}_3  [(d_1 x r_k) matrix]
%                   Xr.res_range = r_{k+1}     [(m x 1) vector].
%                   Xr.U    = U^{(k)}          [(d_1 x r_k) matrix w/
%                                                   orthonormal columns]
%                   Xr.V    = V^{(k)}          [(d_2 x r_k) matrix w/
%                                                   orthonormal columns]
%                   Dense matrix X^{(N)} can be recovered from Xr
%                   using the function 'get_full_mat' such that
%                   X^{(N)} = U^{(k)}*(\tilde{Gamma}_1*V^{(k)'}+\tilde{Gamma}_2)
%                   + \tilde{Gamma}_3*V^{(k)'} + P_{\Omega}^*(r_{k+1}).
% outs: struct. Contains algorithmic information such as the smoothing
%                   parameters eps, level of detail depends on 'saveiterates' 
%                   and 'verbose'. Fields include:
%       - X         (1 x N) cell, with structs: As Xr, but for each
%                   iteration (only if opts.saveiterates == 1).
%       - N         int. Number of performed outer iterations. 
%       - eps       (1 x N) vector. Smoothing parameters '\epsilon_k' for 
%                   each iteration.
%       - r_greatereps (1 x N) vector. Current active ranks 'r_k' for each
%                   iteration, such that r_k = |\{i \in [d]:
%                   \sigma_i(X^{(k)}) > \epsilon_k\}|, cf. [1], [2].
%       - time      (1 x N) vector. Cumulative timestamps at each
%                   iteration.
% =========================================================================
% Further reference (for 'mode_eps'):
% [3] S. Voronin, I. Daubechies, "An iteratively reweighted least squares
% algorithm for sparse regularization", arXiv:1511.08970v3.
% =========================================================================
% Author: Christian Kuemmerle, Johns Hopkins University, kuemmerle@jhu.edu,
% 2019-2020.
% =========================================================================
%% Obtain algorithmic options
[r_upper,N0,N0_inner,p,tol_CG_fac,tol_outer,objective,...
    type_mean,qmean_para,mode_linsolve,mode_eps,epsmin,use_min,...
    increase_antisymmetricweights,tracking,verbose,recsys,saveiterates]=...
    get_options(opts);
%% Obtain problem information
[n,meas_op,y]=get_problem_parameters(prob);
d=min(n.d1,n.d2);
D=max(n.d1,n.d2);
%% Obtain oracle information, if applicable
r_upper = min(r_upper,d);
if not(isfield(prob,'r')) || isempty(prob.r) || not(strcmp(mode_eps,'oracle_model_order'))
    oracle_knowledge = 0;
    r = [];
else
    oracle_knowledge = 1;
    r           = prob.r;        % Rank to be used in the definition of 
                                 % smoothing parameter choice rule
end
%% Initialize logging information
eps             = zeros(1,N0);
r_greatereps    = zeros(1,N0);
time            = zeros(1,N0);
if saveiterates
    X           = cell(1,N0);
    sings       = cell(1,N0);
end
if verbose > 1
    N_inner = zeros(1,N0);
    resvec  = cell(1,N0);
    quot_r  = zeros(1,N0);
    cond_Xell = zeros(1,N0);
    mode_linsolve_it = cell(1,N0);
end 
% =========================================================================
% (only relevant for recommender system application)
if recsys
    clip_ratingscale = opts.clip_ratingscale;
    if clip_ratingscale
        clip_ratingscale = [min(prob.train.values), max(prob.train.values)];
    end
    RMSE         = zeros(1,N0);
    RMSE_train   = zeros(1,N0);
    MABS         = zeros(1,N0);
    MABS_train   = zeros(1,N0);
    prediction.test.row = prob.test.rowind;
    prediction.test.col = prob.test.colind;
    prediction.train.row = prob.train.rowind;
    prediction.train.col = prob.train.colind;
    if verbose > 0
        prediction.test.val = cell(1,N0);
        prediction.train.val = cell(1,N0);
    end
end
% =========================================================================
if tracking  % used if objective value F_{\epsilon,k} is to be tracked
    obj_after_eps_update = zeros(1,N0);
    obj_before_eps_update = zeros(1,N0);
end
%% Determine current rank r_c
if strcmp(objective,'JMLRpaper')
    r_c = d-1;
else
    if oracle_knowledge
        r_c = r;
    else
        r_c = r_upper;
    end
end
%% Initialization
weight_op.U = eye(n.d1,r_c);
weight_op.V = eye(n.d2,r_c);
weight_op.p = p;
eps_c=1;
[weight_op.Wl_inv,weight_op.Wl_inv_eps] = set_weight_infovec(zeros(r_c,1),...
    eps_c,p,objective);
eps_c = Inf;
weight_op.eps = eps_c;
weight_op.sing = zeros(r_c,1);
weight_op_previous = weight_op;
z_c = [];
%% Start iterations
if verbose > 0
    fprintf(['Run MatrixIRLS with p=%.2f, data fit parameter lambda=%.2e and smoothing parameter rule "%-11s"...\n'], ...
    p,lambda,mode_eps);
end
k=1;
first_iterate = 1;
tic
while k <= N0
    %% Calculate the solution of the minimization problem
    if mod(k,25) == 0 
        disp(['Begin Iteration ',num2str(k),'...']);
    end
    eps_previous = eps_c;
    if first_iterate
        X_c_previous.Gam1 = zeros(r_c,r_c);
        X_c_previous.Gam2 = zeros(r_c,n.d2);
        X_c_previous.Gam3 = zeros(n.d1,r_c);
        X_c_previous.res_range   = zeros(1,length(meas_op.Omega));
        X_c_previous.U = weight_op.U;
        X_c_previous.V = weight_op.V;
        X_prev_newbase = zeros(r_c*(r_c+n.d1+n.d2),1);
    else
        X_c_previous   = X_c;
        X_prev_newbase = X_c_newbase;
    end
    %% Update the iterate: Determine X^{(k)}
    [X_c,X_Tn_c,z_c,N_inner_c,resvec_c,relres_c]=update_iterate(weight_op,...
        meas_op,n,y,X_prev_newbase,z_c,...
        lambda,tol_outer,tol_CG_fac,mode_linsolve,increase_antisymmetricweights,...
        first_iterate,N0_inner,...
        type_mean,qmean_para,verbose,objective,X_c_previous.res_range.');
    %% Calculate norm difference of consecutive iterates for stopping criterion (and epsilon rule)
    diff_c  = get_frob_diff_compact(X_c,X_c_previous,meas_op.sps_plc);
    if diff_c < 1e-7 % In this case, the outcome of get_frob_error is not very reliable, therefore compute it again
        Xc_full = get_densemat_from_compact(X_c,meas_op.sps_plc);
        Xc_prev_full = get_densemat_from_compact(X_c_previous,...
            meas_op.sps_plc);
        diff_c = norm(Xc_full-Xc_prev_full,'fro');
    end     
    rel_chg = diff_c/norm_frob_compact(X_c,meas_op.sps_plc);
    if not(first_iterate)
        if rel_chg < tol_outer
            N0 = k;
        end
        if N_inner_c == 0
            N0 = k-1;
        end
    end
    weight_op_previous = weight_op;
    %% Define handle for obtaining spectral information
    X_c_handle = get_matrix_handle(X_c,meas_op,n,objective);
    %% Obtain spectral information
    if strcmp(objective,'pluseps_squared')
        [weight_op.U, sing_c, weight_op.V] = svd(X_c_handle);
        sing_c_mat = zeros(D,D);
        sing_c_mat(1:size(weight_op.U,1),1:size(weight_op.V,1))= sing_c;
        sing_c = sing_c_mat;
    else
        if oracle_knowledge
            [weight_op.U, sing_c, weight_op.V]=bksvd_mod(X_c_handle,...
                min(d,max(r_c,r+1)), 20);
        else
            [weight_op.U, sing_c, weight_op.V]=bksvd_mod(X_c_handle,...
                min(d,r_c), 20);
        end
    end
    weight_op.sing=diag(sing_c);
    %% Initalize current rank
    if first_iterate
        first_iterate = 0;
        eps0 = diff_c;
        if strcmp(objective,'pluseps_squared')
            r_c = max(d1,d2);
        else
            if oracle_knowledge
                r_c = r+1;
            else
                r_c = min(10,r_upper);
            end
        end
    end
    %% Update smoothing parameter \epsilon_k
    eps_c = update_epsilon(eps_c,mode_eps,...
    oracle_knowledge,n,epsmin,use_min,p,weight_op.sing,r,...
    diff_c,eps0,tol_outer,k,N0);
    %% Update weight operator: W^{(k)}
    weight_op.eps = eps_c;
    if not(strcmp(objective,'pluseps_squared'))
        [r_c,weight_op.sing,weight_op.U,weight_op.V] = update_rankpara(X_c_handle,...
            weight_op.sing,weight_op.U,weight_op.V,...
            eps_c,r_upper,r_c,0);
    end
    %%% Calculate quotient quot_r_c between smallest singular value larger 
    %%% than \epsilon_k and \epsilon_k
    if r_c == 0
        N0 = k-1;
        quot_r_c = 1;
    else
        quot_r_c = weight_op.sing(r_c)./eps_c;
    end    
    %% Change of basis for current iterate
    x_c_newbase_resi=proj_tangspace_from_Oprange(X_c.res_range.','MatrixCompletion',...
        weight_op.U,weight_op.V,meas_op.sps_plc,'tangspace');%% update this for increase_antisymmetricweights
    x_c_newbase_Tm1 =proj_tangspace(X_Tn_c,'MatrixCompletion',...
        weight_op_previous.U,weight_op_previous.V,weight_op.U,weight_op.V);
    X_c_newbase=x_c_newbase_resi+x_c_newbase_Tm1;
    
    [weight_op.Wl_inv,weight_op.Wl_inv_eps] = ...
        set_weight_infovec(weight_op.sing(1:r_c),eps_c,p,objective);
    
    if tracking
        %% Track value of objective function
        if strcmp(objective,'pluseps_squared')
            sing_all=weight_op.sing;
        else
            X_c_full = get_densemat_from_iterate(X_c,meas_op.sps_plc);
            [~,singval_mat_c,~]=svd(X_c_full);
            sing_all = diag(singval_mat_c);
        end
        obj_before_eps_update(k) = get_rankobjective(sing_all,eps_previous,p,objective);
        obj_after_eps_update(k) = get_rankobjective(sing_all,eps_c,p,objective);
    end
    
    %% Save output information
    if eps_c <= epsmin
        N0 = k;
    end
    eps(k)=eps_c;
    r_greatereps(k)=r_c;
    time(k) = toc;
    if saveiterates
        %% Save all iterates
        X{k}    = X_c;
        sings{k}= weight_op_previous.sing;
    end
    if verbose > 1
        N_inner(k) = N_inner_c;
        resvec{k}  = resvec_c;
        quot_r(k)  = quot_r_c;
        if r_c == 0
            cond_Xell(k) = 0;
        else
            cond_Xell(k) = weight_op.sing(1)./eps_c;
        end
        mode_linsolve_it{k} = mode_linsolve;
    end
    if recsys 
        %% Calculate train/test error for recommender systems
        if not(isfield(prob,'train') && isfield(prob,'test'))
            error('Please provide training and test data (recommender system application).')
        end
        training_flag = 1;
        prediction_train_val_c = get_prediction_RecSys_IRLS(X_c,...
           prob.train.rowind, prob.train.colind, prob.centscale, clip_ratingscale,training_flag); %prob.train.data
        training_flag = 0;
        prediction_test_val_c = get_prediction_RecSys_IRLS(X_c,...
           prob.test.rowind, prob.test.colind, prob.centscale, clip_ratingscale,training_flag); %prob.test.data
        if verbose > 0
            prediction.train.val{k} = prediction_train_val_c;
            prediction.test.val{k} = prediction_test_val_c;
        end
        RMSE(k) = get_RMSE(prediction_test_val_c,prob.test.values);
        RMSE_train(k) = get_RMSE(prediction_train_val_c,prob.train.values);
        MABS(k) = get_MABS(prediction_test_val_c,prob.test.values);
        MABS_train(k) = get_MABS(prediction_train_val_c,prob.train.values);
        if verbose > 0
            %%% Print algorithmic information
            if k == 1
                fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
                'relres: %.3e  eps_k: %.3e r_k: %d qt_rk: %.2f relchg: %.2e RMSE: %.3f RMSE_train: %.3f\n'], ...
                'init',k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg,RMSE(k),RMSE_train(k)); 
            else
                fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
                'relres: %.3e  eps_k: %.3e r_k: %d  qt_rk: %.2f relchg: %.2e RMSE: %.3f RMSE_train: %.3f\n'], ...
                mode_linsolve,k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg,RMSE(k),RMSE_train(k));  
            end
        end
    else
        if verbose > 0
            %%% Print algorithmic information
            if k == 1
                fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
                'relres: %.3e  eps_k: %.3e r_k: %d qt_rk: %.2f relchg: %.2e\n'], ...
                'init',k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg); 
            else
                fprintf(['%-11s k: %3d  N_inner: %5d   ', ...
                'relres: %.3e  eps_k: %.3e r_k: %d qt_rk: %.2f relchg: %.2e\n'], ...
                mode_linsolve,k,N_inner_c,relres_c,eps_c,r_c,quot_r_c,rel_chg);  
            end
        end
    end
    
    k=k+1;
end
%% Tidy up the size of the relevant output arrays and cells
Xr      = X_c;
if saveiterates
    outs.X      = X(1:N0);
    outs.sings  = sings(1,1:N0);
else
    outs.sings  = weight_op_previous.sing;
end
outs.N      = N0;
outs.eps    = eps(1:N0);
outs.r_greatereps   = r_greatereps(1:N0);
outs.time   = time(1:N0);
if recsys
    outs.RMSE      = RMSE(1:N0);
    outs.RMSE_train= RMSE_train(1:N0);
    outs.MABS      = MABS(1:N0);
    outs.MABS_train= MABS_train(1:N0);
    if verbose > 0 
        prediction.train.val = prediction.train.val(1:N0);
        prediction.test.val = prediction.test.val(1:N0);
        outs.prediction = prediction;
    else
        outs.prediction.train.val = prediction_train_val_c;
        outs.prediction.test.val = prediction_test_val_c;
    end
end
if verbose > 0
    outs.Xr_handle = X_c_handle;
    if verbose > 1
        outs.N_inner  = N_inner(1:N0);
        outs.resvec   = resvec(1:N0);
        outs.quot_r   = quot_r(1:N0);
        outs.cond_Xell = cond_Xell(1:N0);
        outs.mode_linsolve_it = mode_linsolve_it(1:N0);
    end
end
if tracking
    outs.obj_epsA = obj_before_eps_update;
    outs.obj_epsB = obj_after_eps_update;
end
end

function gam_new = handle_M(gam,U,V,weight_vec_c,type_mean,PTPhstarPhPTstar_handle)
d1=size(U,1);
d2=size(V,2);
gam=proj_tangspace(gam,'MatrixCompletion',U,V,U,V);
gam_new=proj_tangspace(PTPhstarPhPTstar_handle(gam) + ...
        apply_weight_vec(gam,d1,d2,weight_vec_c,type_mean),'MatrixCompletion',U,V,U,V);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% To be used in mode == 'tangspace'


function [X_c,X_Tn_c,z_c,N_inner_c,resvec_c,relres_c] = update_iterate(weight_op,...
        meas_op,n,y,...
        x_prev_newbase,z_start,...
        lambda,tol_outer,tol_CG_fac,mode_linsolve,increase_antisymmetricweights,...
        first_iterate,N0_inner,...
        type_mean,qmean_para,verbose,varargin)

U = weight_op.U;
V = weight_op.V;
X_c = struct;
X_c.U = U;
X_c.V = V;
r_c = size(U,2);
d1 = n.d1;
d2 = n.d2;
Wl_inv      = weight_op.Wl_inv;
Wl_inv_eps  = weight_op.Wl_inv_eps;
eps_c       = weight_op.eps;

%%% Update tolerance parameter of the conjugate gradient method.
tol_CG=max(max(min(eps_c.*tol_CG_fac,1e-2)),1e-15);

%%% Load additional variables
if increase_antisymmetricweights
    sing = weight_op.sing;
    p    = weight_op.p;
    if nargin >= 15
        objective = varargin{1};
        X_c4 = varargin{2};
    else
        error('Error in the argument of update_iterate: Please specify "objective" if antisymmetric weights are increased.')
    end
else
    X_c4 = varargin{2};
end

if first_iterate
    X_c.Gam1 = zeros(r_c,r_c);
    X_c.Gam2 = zeros(r_c,d2);
    X_c.Gam3 = zeros(d1,r_c);
    X_c.res_range = Wl_inv_eps.*y./(Wl_inv_eps);
    X_Tn_c = zeros(r_c*(r_c+d1+d2),1);
    z_c = [];
    N_inner_c = 0;
    relres_c = 0;
    resvec_c = 0;
else
    if strcmp(mode_linsolve,'tangspace')
        if increase_antisymmetricweights
            [H,~,~,dH] = weight_op_lowrank_prec(U,V,[Wl_inv;Wl_inv_eps],Wl_inv_eps,1,...
            type_mean,0,[],[],qmean_para,p,sing,eps_c,objective);
        else
            [H,~,~,dH] = weight_op_lowrank_prec(U,V,[Wl_inv;Wl_inv_eps],Wl_inv_eps,1,...
            type_mean,0,[],[],qmean_para);
        end
        PTPhstarPhPTstar_handle = @(gam) PTPhstarPhPTstar_tangspace(gam,U,...
            V,meas_op.sps_plc,meas_op.rowind,meas_op.colind,0);
        if first_iterate
             H{1}=lambda.*H{1};
             dH{1}=lambda.*dH{1}(1:r_c);
             dH{2}=lambda.*dH{2}(1:r_c);
        else
            H{1}=((Wl_inv_eps+lambda)/Wl_inv_eps).*H{1}./(1-H{1});
            dH{1}=((Wl_inv_eps+lambda)/Wl_inv_eps).*dH{1}(1:r_c)./(1-dH{1}(1:r_c));
            dH{2}=((Wl_inv_eps+lambda)/Wl_inv_eps).*dH{2}(1:r_c)./(1-dH{2}(1:r_c));
        end
        weight_vec_c = get_weight_vec(d1,d2,H,dH,increase_antisymmetricweights);
        %%% Define function handle M(gam) = PTPhstarPhPTstar_handle(gam) +  (weight_vec_c.*gam);
        M = @(gam) PTPhstarPhPTstar_handle(gam) + ...
            apply_weight_vec(gam,d1,d2,weight_vec_c,...
            increase_antisymmetricweights,U,V);
        rhs = proj_tangspace_from_Oprange(y.','MatrixCompletion',...
            U,V,meas_op.sps_plc,'tangspace',0);
        prec = @(z) apply_weight_vec(z,d1,d2,(ones(r_c*(d1+d2+r_c),1)+(weight_vec_c)).^(-1),...
            increase_antisymmetricweights);
%         prec = @(z) z;
        if first_iterate
            gamma_start = zeros(r_c*(d1+d2+r_c),1);
        else
            gamma_start = apply_weight_vec(x_prev_newbase,d1,d2,...
                (((Wl_inv_eps+lambda)/Wl_inv_eps)^(-1).*weight_vec_c+ones(r_c*(r_c+d1+d2),1)).^(-1),...
                increase_antisymmetricweights,U,V);
        end
        [gamma,~,relres_c,N_inner_c,resvec_c] = pcg_1(M,rhs,tol_CG,N0_inner,prec,[],gamma_start);
        if first_iterate
            X_Tn_c= gamma;
        else
            X_Tn_c= apply_weight_vec(gamma,d1,d2,...
                (Wl_inv_eps.*weight_vec_c)./(Wl_inv_eps+lambda)+ones(r_c*(r_c+d1+d2),1),...
                increase_antisymmetricweights,U,V);
        end
        res_c =y.'-proj_Oprange_tangspace(gamma,'MatrixCompletion',...
            U,V,meas_op.rowind,meas_op.colind,'tangspace',0);%increase_antisymmetricweights);
        if first_iterate
            X_Tn_c_delta=Wl_inv_eps.*proj_tangspace_from_Oprange(res_c,'MatrixCompletion',...
                U,V,meas_op.sps_plc,'tangspace',0)...;%increase_antisymmetricweights)...
                ./lambda; 
        else
            X_Tn_c_delta=Wl_inv_eps.*proj_tangspace_from_Oprange(res_c,'MatrixCompletion',...
                U,V,meas_op.sps_plc,'tangspace',0)...;%increase_antisymmetricweights)...
                ./(Wl_inv_eps+lambda); 
        end
        X_Tn_c=X_Tn_c-X_Tn_c_delta;

        [X_Tn_c_1,X_Tn_c_2,X_Tn_c_3] = get_Tk_matrices(X_Tn_c,d1,d2,r_c);
        X_c.Gam1 = X_Tn_c_1;
        X_c.Gam2 = X_Tn_c_3;
        X_c.Gam3 = X_Tn_c_2;
        if first_iterate
            X_c.res_range = Wl_inv_eps.*res_c.'./lambda;
        else
            X_c.res_range = Wl_inv_eps.*res_c.'./(Wl_inv_eps+lambda);
        end
        z_c= [];
    else
        error("Other linear system formulations than 'tangspace' currently not available.")
    end
end
end

function eps_c = update_epsilon(eps_c,mode_eps,...
    oracle_knowledge,n,epsmin,use_min,p,varargin)

sing_c   = varargin{1};

if strcmp(mode_eps,'oracle_model_order')
    if oracle_knowledge
        r = varargin{2};
        if p == 1
            eps_rule = sing_c(r+1)/min(n.d1,n.d2);
        else
            eps_rule = sing_c(r+1)/min(n.d1,n.d2)^(p/(2-p));
        end
    else
        error('This epsilon rule is only possible if oracle knowledge is available.')
    end
elseif strcmp(mode_eps,'auto_decay')
        eps0 = varargin{4};
        tol_outer = varargin{5};
        k    = varargin{6};
        N0   = varargin{7};
        
        eps_rule = eps0.*(tol_outer/eps0)^(((k-1)/(N0-1))^(2-p));
elseif strcmp(mode_eps,'iter_diff')
        diff_c = varargin{3};
        
        eps_rule = 2*diff_c/sqrt(min(n.d1,n.d2));
end

if use_min
   eps_c=max(min(eps_rule,eps_c),epsmin); 
else
   eps_c=max(eps_rule,epsmin);
end
end


function X_c_handle = get_matrix_handle(X_c,meas_op,n,objective)
setSval(meas_op.sps_plc,X_c.res_range,length(X_c.res_range));
if strcmp(objective,'pluseps_squared')
    X_c_handle = X_c.U*(X_c.Gam1*weight_op.V'+X_c.Gam2)+X_c.Gam3*X_c.V'+...
        meas_op.sps_plc;
else
    Xc_hdl  = @(x) (X_c.U*(X_c.Gam1*(X_c.V'*x) +(X_c.Gam2*x))...
        + X_c.Gam3*(X_c.V'*x) +meas_op.sps_plc*x);
    Xct_hdl = @(y) (X_c.V*(X_c.Gam1'*(X_c.U'*y)+X_c.Gam3'*y)...
        + X_c.Gam2'*(X_c.U'*y) +meas_op.sps_plc'*y);
    X_c_handle={Xc_hdl,Xct_hdl,n.d1,n.d2};
end
end

function gam_new = PTPhstarPhPTstar_tangspace(gam,U,V,...
            sps_plc,rowind,colind,increase_antisymmetricweights)
m=length(rowind); 
[d1,R] = size(U);
d2= size(V,1);

if increase_antisymmetricweights
    tmp = proj_Oprange_tangspace(gam,'MatrixCompletion',U,V,rowind,colind,...
        'tangspace',increase_antisymmetricweights);
    gam_new = proj_tangspace_from_Oprange(tmp,'MatrixCompletion',...
                    U,V,sps_plc,'tangspace',increase_antisymmetricweights);
else
    tmp = proj_Oprange_tangspace(gam,'MatrixCompletion',U,V,rowind,colind,'tangspace');
    
    gam_new = proj_tangspace_from_Oprange(tmp,'MatrixCompletion',...
                    U,V,sps_plc,'tangspace');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_upper,N0,N0_inner,p,tol_CG_fac,tol_outer,objective,...
    type_mean,qmean_para,mode_linsolve,mode_eps,epsmin,use_min,...
    increase_antisymmetricweights,tracking,verbose,recsys,saveiterates]=...
    get_options(opts)

if isfield(opts,'R')
    r_upper = opts.R; % Upper bound for number of singular values to be calculated
else
    r_upper = Inf;
end


N0     = opts.N0;       % Maximal number of iterations for the IRLS algorithm
N0_inner=opts.N0_inner; % Maximal number of inner iterations (e.g. in the CG method)
p      = opts.p;        % Parameter describing the non-convexity of the objective function
                        % (p = 0 for logarithmic objective, 0 < p < 1 for
                        % Schatten-p quasi-norm objective)
tol_CG_fac    = opts.tol_CG_fac;   % Tolarance used in the stopping criterion of the inner CG method
tol_outer = opts.tol;      % Tolerance used as stopping criterion of outer IRLS iterations
                        % (min. change of relative Frobenius norm)
if not(isfield(opts,'objective')) || isempty(opts.objective)
    objective = 'objective_thesis';
else
    objective   = opts.objective;% objective = 'objective_thesis': Corresponds to using objectives
end
%                           of the type:
%                           log(max(sigma,eps))+0.5.*(min(sigma,eps)^2/eps^2-1).
%                           (as used in the Ph.D. thesis [Kuemmerle
%                           2019])
                        % mode = 'pluseps': Corresponds to using objectives
%                           of the type:
%                           log(max(sigma,eps)+eps)+0.25.*(min(sigma,eps)^2/eps^2-1)
                        % mode = 'pluseps_squared':Corresponds to using objectives
%                           of the type:
%                           log(max(sigma,eps)^2+eps^2)+0.5.*(min(sigma,eps)^2/eps^2-1)
type_mean = opts.type_mean;% type_mean = 'harmonic': corresponds to usage 
                        % of harmonic mean of weight matrices, as Kuemmerle
if strcmp(type_mean,'qmean')
    qmean_para = opts.qmean_para; 
    % parameter for the q-mean: = 1 for arithmetic mean,
    %                           = 0 for heometric mean,
    %                           = -1 for harmonic mean,
    %                           = -infty for "min" mean.
elseif strcmp(type_mean,'optimal')
    if p == 0
        type_mean = 'geometric';
        qmean_para = [];
    else
        type_mean = 'qmean';
        qmean_para = p/(p-2);
    end
else
    qmean_para = [];
end
mode_linsolve = opts.mode_linsolve;

 %(just incorporated for mode = 'fullspace' yet)    % and Sigl (2017).
                        %    type_mean = 'Wasserstein': corresponds to a
                        %    Wasserstein mean of the inverses of the
                        %    one-sided weight matrix inverses, cf. Bhatia,
                        %    Jain, Lim (2018)
mode_eps   = opts.mode_eps;%= 'classical_r': eps choice dependent on true
%                                           rank r, as classical
%                          = 'auto_decay' : eps that is automatically
%                             superlinearly decaying, rule such that eps0=
%                             sing_val_1(x^(1)) and 
%                             eps(iter) = 10^(-(iter/N0)^(2-p)*10)
%                          = 'iter_diff'  : eps similar to Voronin &
%                                           Daubechies (2016), i.e.,
%                             eps(iter)= min(eps(iter-1),\|x{iter}-x{iter-1}\|_2
%                                      + eps0*10^(-(iter/10)*(2-p)*10))
epsmin     = opts.epsmin;  
use_min    = opts.use_min; % == 1: In the epsilon choice, use the minimum 
                           %        of the previous one and the one calculated
                           %        by the rule.
                           % == 0: Use the one calculated by the epsilon rule.
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    saveiterates = 1;
else
    saveiterates = 0;
end

if isfield(opts,'recsys') && opts.recsys == 1
    recsys = 1;
else
    recsys = 0;
end

increase_antisymmetricweights = opts.increase_antisymmetricweights;


if isfield(opts,'tracking') &&  ~isempty(opts.tracking)
    tracking=opts.tracking;
else
    tracking = 0;
end

verbose     = opts.verbose; % verbose = 0: No output except of last iterate
                            % verbose = 1: Return also intermediate
                            % iterates and a few algorithmic parameters 
                            % verbose = 2: Return also many further
                            % parameters

end

function  [n,meas_op,Yval]=get_problem_parameters(prob)
[rowind,colind] = find(prob.Phi); 
Yval   = prob.y.';   
Omega  = find(prob.Phi);
n.d1   = prob.d1;       % First dimension of matrix X to be recovered
n.d2   = prob.d2;       % Second dimension of matrix X to be recovered

meas_op.Omega  = Omega;
meas_op.rowind = rowind;
meas_op.colind = colind;
meas_op.sps_plc = sparse(rowind,colind,ones(length(Omega),1),n.d1,n.d2);
end