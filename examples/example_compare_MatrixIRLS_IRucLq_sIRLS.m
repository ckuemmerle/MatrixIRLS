% Example script that compares the iteratively reweighted least squares
% algorithms 'MatrixIRLS' of [1,2] with algorithms 'IRLS-p' and 'sIRLS-p' of
% the paper [3, Mohan and Fazel 2012], and the with the algorithm
% 'tIRucLq'/'IRucLq' of the paper [4, Lai, Xu and Yin 2013], for both the
% non-convexity parameters q (resp. p) = 0 and = 0.5.
%
% For a description of the random model on the sampled entries, see
% documentation of "sample_phi_MatrixCompletion.m". For a description of 
% the model for the matrix X0 to be completed, see documentation of 
% "sample_X0_lowrank.m".
%
% If add_noise == 0, no noise is added to the measurements. As the
% algorithms of [4] are designed to solve an unconstrained problem with
% objective (case q > 0)
% J(X,epsilon) := (1/q).*\sum_{i=1}^{d_1} (\sigma_i^2(X)+epsilon^2)^(q/2) +
%                  1/(2*lambda).*\|P_{\Omega}(X)-y\|_2^2,
% we choose lambda = 1e-6.*norm(y) for 'tIRucLq' and 'IRucLq' in case that
% the prescribed lambda is zero, i.e., if 'opts_custom.lambda == 0'.
% =========================================================================
% Author: Christian Kuemmerle, 2020.
% =========================================================================
% References:
% [1] C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples", ICML 2021.
%
% [2] C. Kuemmerle, C. Mayrink Verdun, "Escaping Saddle Points in 
% Ill-Conditioned Matrix Completion with a Scalable Second Order Method", 
% ICML 2020 Workshop "Beyond First Order Methods in ML Systems".
%
% [3] Karthik Mohan and Maryam Fazel, "Iterative reweighted algorithms for 
% matrix rank minimization", Journal of Machine Learning Research, 
% 13(1):3441-3473, 2012.
%
% [4]  Ming-Jun Lai, Yangyang Xu and Wotao Yin, "Improved iteratively 
% reweighted least squares for unconstrained smoothed \ell_q minimization",
% SIAM Journal on Numerical Analysis 51(2):927-957, 2013.
% =========================================================================
rng('shuffle')
%% Set parameters
% Number of rows of the matrix that defines the problem:
d1 = 400; 
% Number of columns:
d2 = 400; 
% Rank:
r = 15;
% Number of degrees of freedom of the setting:
df_LR = @(rr) rr*(d1 + d2 - rr);
df_LR_val=df_LR(r);
% Oversampling factor:
oversampling = 2.0;
% Number of measurements:
m = floor(min(oversampling*df_LR_val,d1*d2));
%% Sample the measurement matrix Phi (pattern of revealed entries) 
max_nr_resample = 1000;
[Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
[rowind,colind] = ind2sub([d1,d2],Omega);
%% Sample the ground truth matrix X0 and measured entries y
modeX0      = 'condition_control_log'; % 1 , 2.
complexflag = 0;
%%%%%%%%%%%%%%%%%%%%%%
%%%% Set condition number:
cond_nr = 5e0;
%%%%%%%%%%%%%%%%%%%%%%
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
X0 = {U0,V0};
% X0_full=U0*V0';
y = partXY(U0',V0',rowind,colind,m).';
%% If desired: add noise 
add_noise = 0;
noise_fac = 0.1; % Multiplicative factor determining noise level of additive
               % noise
if add_noise
    noise = randn(m,1); 
    yn = y + noise_fac./norm(y).*noise./norm(noise,2);
else
    yn = y;
end

%% Choose algorithms for matrix completion
alg_names={'MatrixIRLS','tIRucLq','sIRLS-p','IRLS-p'}; %'IRucLq'
opts_custom.p = [0]; % (non-)convexity parameters p for IRLS
% p = 0: sum of log objective
% 0 < p < 1: Schatten-p quasi-norm.
% p = 1: Nuclear norm.
%%% Set optional algorithmic parameters
opts_custom.N0 = 400;           % Max. number of (outer) iterations.
opts_custom.N0_inner = 200;     % Max. number of (inner) iterations for MatrixIRLS.
opts_custom.tol_CG_fac=1e-4;    % tolerance for stopping criterion of inner iterations
opts_custom.type_mean = {'optimal'};        % 'harmonic','arithmetic','geometric','min'
opts_custom.objective = 'objective_thesis'; % 'pluseps','pluseps_squared'
opts_custom.mode_eps = 'oracle_model_order';%'oracle_model_order';% 'iter_diff','auto_decay'
opts_custom.lambda = 0;
opts_custom.saveiterates = 1;
opts_custom.verbose = 1;

%% Run algorithms for matrix completion
[Xr,outs,alg_names] = run_MC_algos(Phi,yn,r,alg_names,opts_custom);

%% 
plotflag =1;
[sings,singsX0] = get_singvals_algorec(Xr,X0,Phi,r,'last',...
plotflag,alg_names);
%% Calculating the error statistics
verbose = 1; % provide text output after error calculations (or not if = 0)
frob_mode = 'full';
[error_fro_rel,~] = get_frob_errors(Xr,X0,Phi,alg_names,...
    frob_mode,verbose);

%% Visualization of Frobenius errors of iterates of algorithms
visualize = 1;
if visualize == 1
    visualize_errorcurves_combined(error_fro_rel,alg_names);
    plot_times_errors(outs,error_fro_rel,alg_names);
end


%% Run everything again for p = 0.5
alg_names={'MatrixIRLS','tIRucLq','sIRLS-p','IRLS-p'};
opts_custom.p = [0.5];
Phi = sparse(rowind,colind,ones(m,1),d1,d2);
%% Run algorithms for matrix completion
[Xr,outs,alg_names] = run_MC_algos(Phi,yn,r,alg_names,opts_custom);
plotflag =1;
[sings,singsX0] = get_singvals_algorec(Xr,X0,Phi,r,'last',...
plotflag,alg_names);
frob_mode = 'full';
[error_fro_rel,~] = get_frob_errors(Xr,X0,Phi,alg_names,...
    frob_mode,verbose);
if visualize == 1
    visualize_errorcurves_combined(error_fro_rel,alg_names);
    plot_times_errors(outs,error_fro_rel,alg_names);
end