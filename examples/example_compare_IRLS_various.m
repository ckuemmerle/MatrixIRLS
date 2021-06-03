% Example script that compares the iteratively reweighted least squares
% algorithms 'MatrixIRLS' of [1,2], a similar fast implementation of IRLS
% with suboptimal weight operator 'arithmetic', and (a variant of) 
% the algorithm of [3]  based on row/column reweighting 
% ('IRLS-col','IRLS-row', generalized for 0 <= p <= 1 arbitrary).
%
% For a description of the random model on the sampled entries, see
% documentation of "sample_phi_MatrixCompletion.m". For a description of 
% the model for the matrix X0 to be completed, see documentation of 
% "sample_X0_lowrank.m".
% =========================================================================
% Author: Christian Kuemmerle, 2020-2021.
% =========================================================================
% References:
% [1] C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples", ICML 2021.
%
% [2] C. Kuemmerle, C. Mayrink Verdun, "Escaping Saddle Points in 
% Ill-Conditioned Matrix Completion with a Scalable Second Order Method", 
% ICML 2020 Workshop "Beyond First Order Methods in ML Systems".
%
% [3] M. Fornasier, H. Rauhut, R. Ward, "Low-rank matrix recovery via 
% iteratively reweighted least squares minimization. SIAM J. Optim., 
% 21(4):1614?1640, 2011.

rng('shuffle')
%% Set parameters
% Number of rows of the matrix that defines the problem:
d1 = 400; 
% Number of columns:
d2 = 400; 
% Rank:
r = 5;
% Number of degrees of freedom of the setting:
df_LR = @(rr) rr*(d1 + d2 - rr);
df_LR_val=df_LR(r);
% Oversampling factor:
oversampling = 4.0;
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
cond_nr = 2e0;
%%%%%%%%%%%%%%%%%%%%%%
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
X0 = {U0,V0};
% X0_full=U0*V0';
y = partXY(U0',V0',rowind,colind,m).';
%% If desired: add noise 
add_noise = 0;
noise_fac = 1; % Multiplicative factor determining noise level of additive
               % noise
if add_noise
    noise = randn(m,1); 
    yn = y + noise_fac.*noise./norm(noise,2);
else
    yn = y;
end

%% Choose algorithms for matrix completion
alg_names={'MatrixIRLS','IRLS-row','IRLS-col'};

%%% Set optional algorithmic parameters
opts_custom.tol = 1e-7;        % tolerance for stopping criterion
opts_custom.N0 = 400;           % Max. number of (outer) iterations for 
                                % 'second-order algorithms', which include 
                                % MatrixIRLS, R2RILS and RTRMC.
opts_custom.N0_firstorder = 1000; % Max. number of iterations for 'first-order algorithms'.
%%% Optional parameters addressing options of 'second-order algorithms'.
opts_custom.tol_CG = 1e-4;    % tolerance for stopping criterion of inner iterations
opts_custom.N0_inner = 500;     % Max. number of (inner) iterations for 'second-order algorithms'


%%% Optional parameters addressing only options for IRLS
opts_custom.p = [0]; % (non-)convexity parameters p for IRLS
% p = 0: sum of log objective
% 0 < p < 1: Schatten-p quasi-norm.
% p = 1: Nuclear norm.
opts_custom.R = min(d1,d2);             %min(d1,d2);%floor(10*r);50;%
opts_custom.adaptive_cg = 0;
opts_custom.mode_linsolve = 'tangspace';    % 'rangespace'
opts_custom.type_mean = {'optimal','arithmetic'};        % 'harmonic','arithmetic','geometric','min'
opts_custom.objective = 'objective_thesis'; % 'pluseps','pluseps_squared'
opts_custom.mode_eps = 'oracle_model_order';% 'iter_diff','auto_decay'
opts_custom.epsmin = 1e-16;
opts_custom.tracking = 0;
opts_custom.lambda = 0;
opts_custom.increase_antisymmetricweights=0;
opts_custom.saveiterates = 1;
opts_custom.verbose = 1;

%% Run algorithms for matrix completion
[Xr,outs,alg_names] = run_MC_algos(Phi,yn,r,alg_names,opts_custom);

%% Calculating the error statistics
calculate_partial_Frob_norms = 0; % if =1, calculate also partial Frobenius norms (restricted to Omega / Omega^c)
verbose = 1; % provide text output after error calculations (or not if = 0)
frob_mode = 'full';
[error_fro_rel,error_fro] = get_frob_errors(Xr,X0,Phi,alg_names,...
    frob_mode,verbose);
if calculate_partial_Frob_norms
    frob_mode = 'Phi';
    [erPhi_fro_rel,erPhi_fro] = get_frob_errors(Xr,X0,Phi,alg_names,...
        frob_mode,verbose);
    frob_mode = 'Phi_comp';
    [erPhic_fro_rel,erPhic_fro] = get_frob_errors(Xr,X0,Phi,alg_names,...
        frob_mode,verbose);
end

%% Visualization of Frobenius errors of iterates of algorithms
visualize = 1;
if visualize == 1
    visualize_errorcurves_combined(error_fro_rel,alg_names);
    plot_times_errors(outs,error_fro_rel,alg_names);
end