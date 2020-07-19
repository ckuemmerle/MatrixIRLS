% This script corresponds to the Figure 5 contained in the paper "Escaping 
% Saddle Points in Ill-Conditioned Matrix Completion with a Scalable
% Second Order Method" by Christian Kuemmerle and Claudio M. Verdun,
% published at the ICML 2020 Workshop on 'Beyond First Order Methods in 
% Machine Learning'.

% It performs a completion task for a highly ill-conditioned 1000 x 1000
% matrix of rank r=10 with kappa=1e10 and with an oversampling rate of
% rho=4.
%
% Author: Claudio M. Verdun, June 2020.
%% Set parameters
% Number of rows of the matrix that defines the problem:
d1 = 1000; 
% Number of columns:
d2 = 1000; 
% Rank:
r = 10;
% Number of degrees of freedom of the setting:
df_LR = @(rr) rr*(d1 + d2 - rr);
df_LR_val=df_LR(r);
% Number of measurements:
m = floor(min(4*df_LR_val,d1*d2));
rng('shuffle')
%% Sample the measurement matrix Phi (pattern of revealed entries) 
max_nr_resample = 1000;
[Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
[rowind,colind] = ind2sub([d1,d2],Omega);
%% Sample the ground truth matrix X0 and measured entries y
modeX0      = 'condition_control_log'; % Singular values are interpolated 
% in a logarithmic way between 1 and kappa
complexflag = 0;
cond_nr = 1e10; % condition number of the low-rank matrix X0
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
X0 = {U0,V0};
y = partXY(U0',V0',rowind,colind,m).';
%% Choose algorithms for matrix completion
alg_names={'MatrixIRLS','R2RILS','RTRMC','LRGeomCG','LMaFit','ScaledASD','ScaledGD','NIHT'};
opts_custom.p=0; %%% choose (non-)convexity parameter p for IRLS
p=opts_custom.p;
opts_custom.tol = 1e-11;
opts_custom.tol_CG_fac = 1e-5*cond_nr^(-1);
opts_custom.N0 = 400;
opts_custom.N0_firstorder = 4000;
opts_custom.N0_inner = 500; 
opts_custom.mode_linsolve = 'tangspace'; 
opts_custom.type_mean={'geometric'};
opts_custom.epsmin = 1e-16;
opts_custom.lambda = 0;
opts_custom.saveiterates = 1;
opts_custom.verbose = 1;
%% Run algorithms for matrix completion
[Xr,outs] = run_MC_algos(Phi,y,r,alg_names,opts_custom);

%% Calculating the error statistics
verbose = 1; % provide text output after error calculations
frob_mode = 'full';
[error_fro_rel,error_fro] = get_frob_errors(Xr,X0,Phi,alg_names,frob_mode,verbose);

%% Visualization of Frobenius errors vs. times
% visualize_errorcurves_combined(error_fro_rel,alg_names);
plot_options.markers={'-x', '-+', '-*', '-o','-x', '-s', '-d', '-^'};
plot_options.markerssize={5,5,5,5,5,5,5,5};
plot_options.maxtime = 60;
fig = plot_times_errors(outs,error_fro_rel,alg_names,plot_options);
plot_options.ColorOrderIndices=[1,2,3,4,5,6,7,8,6];
postprocess_fig(alg_names,'Time in seconds','Rel. Frob. error of $\textbf{X}^{(k)}$',plot_options);
%% Save the data from simulation
outs = remove_intermediate_iterates(outs);
clear Xr;
filename=strcat('experiment_MCalgos_ICML2020_Fig5_1000x1000_kappa1e10_rank10');
save(strcat(filename,'.mat'))
savefig(strcat(filename,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat(filename,'.tex'))
end