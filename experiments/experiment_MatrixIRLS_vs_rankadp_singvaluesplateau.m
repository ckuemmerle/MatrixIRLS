% This script reproduces the experiment of Figure 5 of the paper
% contained in the paper "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples" by Christian Kuemmerle 
% and Claudio Mayrink Verdun.

% It performs a completion task for a highly ill-conditioned 1000 x 1000
% matrix of rank r=30 with kappa=1e10 from few samples (oversampling factor
% of 1.5 with respect to the degrees of freedom of low-rank matrix).
%
% Authors: Claudio M. Verdun and Christian Kuemmerle, February 2021.

rng('shuffle')
%% Set parameters
% Number of rows of the matrix that defines the problem:
d1 = 1000; 
% Number of columns:
d2 = 1000; 
% Rank:
r = 30;
k = r;%100;
% Number of degrees of freedom of the setting:
df_LR = @(rr) rr*(d1 + d2 - rr);
df_LR_val=df_LR(r);
% Oversampling factor:
oversampling = 1.50;
% Number of measurements:
m = floor(min(oversampling*df_LR_val,d1*d2));
%% Sample the measurement matrix Phi (pattern of revealed entries) 
max_nr_resample = 1000;
[Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
[rowind,colind] = ind2sub([d1,d2],Omega);
%% Sample the ground truth matrix X0 and measured entries y
modeX0      = 'condition_control_log_plateau';
complexflag = 0;
cond_nr = 1e10;
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
X0 = {U0,V0};
y = partXY(U0',V0',rowind,colind,m).';
%% If desired: add noise (optional, is not part of the experiment)
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
alg_names={'LRGeomCG_pursuit','MatrixIRLS','R3MC-rankupd','R3MC','LRGeomCG'};
%%% Set optional algorithmic parameters
opts_custom.tol = 1e-9;             % tolerance for stopping criterion
opts_custom.N0 = 400;               % Max. number of (outer) iterations for 
                                    % 'second-order algorithms', which include
                                    % MatrixIRLS, R2RILS and RTRMC.
                                
opts_custom.N0_firstorder = 4000;   % Max. number of iterations for 'first-order algorithms'.

%%% Optional parameters addressing options of 'second-order algorithms'.
opts_custom.N0_inner = 20;         % Max. number of (inner) iterations for 'second-order algorithms'
%%% Optional parameters addressing only options for IRLS
opts_custom.p = 0; % (non-)convexity parameters p for IRLS
opts_custom.saveiterates = 1;
opts_custom.N_SVD = 3;
opts_custom.tol = 1e-13;
%%% Optional parameter set for R3MC
opts_custom.beta_type = 'P-R'; % Algorithm type. Options: H-S (Hestenes-Stiefel), F-R (Fletcher–Reeves) and P-R (Polak–Ribière)
%% Run algorithms for matrix completion
[Xr,outs,alg_names] = run_MC_algos(Phi,yn,k,alg_names,opts_custom);

%% Calculating the error statistics
verbose = 1;                                        % provide text output after error calculations (or not if = 0)
frob_mode = 'full';
[error_fro_rel,error_fro] = get_frob_errors(Xr,X0,Phi,alg_names,...
    frob_mode,verbose);
%% Visualization of Frobenius errors of iterates of algorithms

curdate = datestr(now,'yyyy-mm-dd_HH-MM-SS');
filename=strcat('MIRLSvsRankadp_',curdate,'_kappa1e10_r30');

visualize = 1;
if visualize == 1
    visualize_errorcurves_combined(error_fro_rel,alg_names);
    plot_times_errors(outs,error_fro_rel,alg_names);
    savefig(strcat('results/',filename,'.fig'))
end
%% Get reconstructions
X = get_densemat_from_iterate(Xr,Phi);
X0_full = X0{1}*X0{2}';
%% Get and visualize singular values
plotflag =1;
[sings,singsX0] = get_singvals_algorec(Xr,X0,Phi,r,'last',...
plotflag,alg_names);
savefig(strcat('results/',filename,'_singvals.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',filename,'.tex'))
end