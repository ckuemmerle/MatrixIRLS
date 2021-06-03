% This script plots the experimental data of Figure 1 contained in the
% paper
% C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples", ICML 2021.
%
% In the experiment, the performance of several algorithms for exact matrix
% completion of (1000 x 1000) matrices of rank 5 with condition number
% kappa = 10 is compared in terms of median rel. Frobenius error of
% reconstruction vs. oversampling factor rho (which is the quotient between
% sample size and number of degrees of freedom).
%
% To re-run the experiments, see file
% "run_experiment_ICML2021_Figure1.m".
clear all;
custom_options = struct;
custom_options.ColorOrderIndices=[2,3,4,5,6,7,8,9,6];

load(['experiment_MCalgos_ICML2020_Fig1_samplesvsmedian_kappa10.mat']);
median_err_rel_experiment1 = median_err_rel;
alg_names_experiment1 = alg_names;
load(['MCsamplecomp_d1_1000_d2_1000_r_5_kappa_10_2021-05-31_20-43-24_rho1.05-4.mat'])
median_err_rel = [median_err_rel_experiment1;median_err_rel];
alg_names = [alg_names_experiment1,alg_names(1)];


plot_results_flexible(median_err_rel,parameters,alg_names,'logy',...
['Median of rel. Frob. errors of $\mathbf{X}^{(K)}$'],...
'Oversampling factor $\rho$',custom_options);
set(gcf,'Position',[744,820,384,230]);

save_filename=strcat('experiment_MCalgos_ICML2021_Fig1_samplesvsmedian_kappa10');
savefig(strcat('results/',save_filename,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',save_filename,'.tex'),...
        'height', '\figureheight', 'width', '\figurewidth',...
        'extraaxisoptions',['xlabel style={font=\tiny},',...
        'ylabel style={font=\tiny},',...
        'legend style={font=\fontsize{7}{30}\selectfont, anchor=south, legend columns = 3, at={(0.5,1.03)}}']);
end