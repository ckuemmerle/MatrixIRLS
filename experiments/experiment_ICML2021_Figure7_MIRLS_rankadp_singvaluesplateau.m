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
alg_names={'LRGeomCG','LRGeomCG_pursuit','R3MC','R3MC-rankupd','MatrixIRLS'}; %'LRGeomCG',
%%% Set optional algorithmic parameters
opts_custom.tol = 1e-13;             % tolerance for stopping criterion
opts_custom.tol_CG = 1e-2;
opts_custom.N0 = 400;
opts_custom.N0_firstorder = 4000;
opts_custom.N0_inner = 500; 
opts_custom.N_SVD = 3;
opts_custom.saveiterates = 1;
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
filename = 'experiment_MCalgos_ICML2021_Fig7_1000x1000_kappa1e10_rank30_rankada';
save(strcat('results/',filename,'.mat'))
% filename=strcat('MIRLSvsRankadp_',curdate,'_kappa1e10_r30');

plot_options = struct;
plot_options.markers={'-o',':>','-v',':*','-x'};%{'-x', '-+', '-*', '-o','-x', '-s', '-d', '-^','-v'};
plot_options.markerssize={7,5,5,3,5};
plot_options.maxtime = 180;
plot_options.ColorOrderIndices=[2,4,3,9,1];

alg_names{2} = 'LRGeomCG Pursuit';
alg_names{4} = 'R3MC w/ Rank Upd';
perm=[1,3,5,2,4];
plot_options.ColorOrderIndices = plot_options.ColorOrderIndices(perm);
plot_options.markers = plot_options.markers(perm);
plot_options.markerssize = plot_options.markerssize(perm);
plot_times_errors(outs(perm),error_fro_rel(perm),alg_names(perm),plot_options);
% plot_options.ColorOrderIndices=[1,2,3,4,5,6,7,9,6];
postprocess_fig(alg_names,'Time in seconds','Rel. Frob. error',plot_options);
% visualize = 1;
% if visualize == 1
%     visualize_errorcurves_combined(error_fro_rel,alg_names);
% plot_times_errors(outs,error_fro_rel,alg_names);
savefig(strcat('results/',filename,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',filename,'.tex'),...
                'height', '\figureheight', 'width', '\figurewidth',...
        'extraaxisoptions',['xlabel style={font=\tiny},',...
        'ylabel style={font=\tiny},',...
        'legend style={font=\fontsize{7}{30}\selectfont, anchor=south, legend columns = 3, at={(0.4,1.03)}}'])
end
% end
%% Get reconstructions
X = get_densemat_from_iterate(Xr,Phi);
X0_full = X0{1}*X0{2}';
%% Visualize singular values
filename_singvals = strcat(filename,'_singvals');
plotflag = 2;
perm=[1,3,5,2,4,6];
plot_options.markers={'o','>','v','*','x','s'};
plot_options.markerssize={150,150,150,140,150,200};
plot_options.ColorOrderIndices=[2,4,3,9,1,8];
plot_options.ColorOrderIndices = plot_options.ColorOrderIndices(perm);
plot_options.markers = plot_options.markers(perm);
plot_options.markerssize = plot_options.markerssize(perm);
alg_names_X0 = [alg_names,'X^0'];
[sings,singsX0] = get_singvals_algorec(Xr(perm(1:5)),X0,Phi,r,'last',...
plotflag,alg_names_X0(perm));
postprocess_fig(alg_names_X0(perm),'Singular value index $i$','Singular value $\sigma_i(X^{(K)})$',plot_options);
savefig(strcat('results/',filename_singvals,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',filename_singvals,'.tex'),...
                'height', '\figureheight', 'width', '\figurewidth',...
...%                 'ymode','log',...
        'extraaxisoptions',['xlabel style={font=\tiny},',...
        'ylabel style={font=\tiny},',...
        'legend style={font=\fontsize{7}{30}\selectfont, anchor=south, legend columns = 3, at={(0.5,1.03)}}'])
end
%% Visualize relative errors on singular values
plotflag = 1;
perm=[1,3,5,2,4];
plot_options.markers={'o','>','v','*','x'};
plot_options.markerssize={150,150,150,140,150};
plot_options.ColorOrderIndices=[2,4,3,9,1];
plot_options.ColorOrderIndices = plot_options.ColorOrderIndices(perm);
plot_options.markers = plot_options.markers(perm);
plot_options.markerssize = plot_options.markerssize(perm);
[~,~] = get_singvals_algorec(Xr(perm),X0,Phi,r,'last',...
plotflag,alg_names(perm));
% plot_options.markers = plot_options.markers(perm);
% plot_options.markerssize = plot_options.markerssize(perm);
postprocess_fig(alg_names(perm),'Singular value index $i$','Rel. error $\frac{|\sigma_i(X^{(K)})-\sigma_i(X^{0})|}{\sigma_i(X^{0})}$',plot_options);
filename_singvalerrors = strcat(filename,'_singvalerrors');
savefig(strcat('results/',filename_singvalerrors,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',filename_singvalerrors,'.tex'),...
                'height', '\figureheight', 'width', '\figurewidth',...
...%                 'ymode','log',...
        'extraaxisoptions',['xlabel style={font=\tiny},',...
        'ylabel style={font=\tiny},',...
        'legend style={font=\fontsize{7}{30}\selectfont, anchor=south, legend columns = 3, at={(0.5,1.03)}}'])
end