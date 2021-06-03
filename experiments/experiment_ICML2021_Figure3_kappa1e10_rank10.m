% This script reproduces the experiment of Figure 4 of the paper
% C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples", ICML 2021.
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
cond_nr = 1e10; % condition number of the low rank matrix
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
X0 = {U0,V0};
y = partXY(U0',V0',rowind,colind,m).';
%% Choose algorithms for matrix completion
alg_names={'MatrixIRLS','R2RILS','RTRMC','LRGeomCG','LMaFit','ScaledASD','ScaledGD','NIHT','R3MC'};
opts_custom.p=0; %%% choose (non-)convexity parameters p for IRLS
p=opts_custom.p;
opts_custom.tol = 1e-12;
opts_custom.tol_CG = 1e-9;
opts_custom.N0 = 400;
opts_custom.N0_firstorder = 4000;
opts_custom.N0_inner = 500; 
opts.N_SVD          = 100;
opts_custom.mode_linsolve = 'tangspace'; 
opts_custom.type_mean={'geometric'};
opts_custom.epsmin = 1e-16;
opts_custom.lambda = 0;
opts_custom.saveiterates = 1;
opts_custom.verbose = 2;
%% Run algorithms for matrix completion
[Xr,outs] = run_MC_algos(Phi,y,r,alg_names,opts_custom);

%% Calculating the error statistics
verbose = 1;
frob_mode = 'full';
[error_fro_rel,error_fro] =get_frob_errors(Xr,X0,Phi,alg_names,frob_mode,verbose);

%% Visualization of Frobenius errors of iterates of algorithms
%     visualize_errorcurves_combined(error_fro_rel,alg_names);
plot_options.markers={'-x', '-+', '-*', '-o','-x', '-s', '-d', '-^','-v'};
plot_options.markerssize={5,5,5,5,5,5,5,5,5};
plot_options.maxtime = 60;
plot_times_errors(outs,error_fro_rel,alg_names,plot_options);
plot_options.ColorOrderIndices=[1,2,3,4,5,6,7,9,6];
postprocess_fig(alg_names,'Time in seconds','Rel. Frob. error',plot_options);
set(gcf,'name','Algorithmic comparisons: Completion of a $1000 \times 1000$ rank-$10$ matrix, $\kappa=10^5$')
%% Save the data from simulation
outs = remove_intermediate_iterates(outs);
clear Xr;
filename=strcat('experiment_MCalgos_ICML2021_Fig3_1000x1000_kappa1e5_rank10');
save(strcat('results/',filename,'.mat'))
savefig(strcat('results/',filename,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',filename,'.tex'),...
                'height', '\figureheight', 'width', '\figurewidth',...
        'extraaxisoptions',['xlabel style={font=\tiny},',...
        'ylabel style={font=\tiny},',...
        'legend style={font=\fontsize{7}{30}\selectfont, anchor=south, legend columns = 3, at={(0.5,1.03)}}']);
end