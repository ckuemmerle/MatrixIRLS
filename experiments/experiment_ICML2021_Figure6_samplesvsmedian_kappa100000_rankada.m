% This script plots the experimental data of Figure 6 contained in the
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
% "run_experiment_ICML2021_Figure6.m".
clear all;
custom_options = struct;
custom_options.ColorOrderIndices=[4,9,2];
custom_options.markers = {'-x', ':>',':*'};



load(['experiment_MCalgos_ICML2020_Fig2_samplesvsmedian_kappa100000.mat']);
median_err_rel_exp1 = median_err_rel([1],:);
alg_names_exp1 = alg_names([1]);
load(['MCsamplecomp_d1_1000_d2_1000_r_5_kappa_100000_2021-05-31_20-57-28_rho1.05-4.mat'])
median_err_rel_combined(1,:) =  median_err_rel_exp1(1,:);
median_err_rel_combined(2,:) = median_err_rel(1,:);
alg_names_combined = alg_names_exp1(1);
alg_names_combined(2) = alg_names(1);
alg_names_combined{2} = 'LRGeomCG Pursuit';
load(['MCsamplecomp_d1_1000_d2_1000_r_5_kappa_100000_2021-05-31_19-57-08_rho1.05-4.mat'])
median_err_rel_combined(3,:) = median_err_rel(2,:);
alg_names_combined(3) = alg_names(2);
alg_names_combined{3} = 'R3MC w/ Rank Upd';

plot_results_flexible(median_err_rel_combined,parameters,alg_names_combined,'logy',...
['Median of rel. Frob. errors of $\mathbf{X}^{(K)}$'],...
'Oversampling factor $\rho$',custom_options);
set(gcf,'Position',[744,820,384,230]);

save_filename=strcat('experiment_MCalgos_ICML2021_Fig6_samplesvsmedian_kappa100000_rankada');
savefig(strcat('results/',save_filename,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',save_filename,'.tex'),...
        'height', '\figureheight', 'width', '\figurewidth',...
        'extraaxisoptions',['xlabel style={font=\tiny},',...
        'ylabel style={font=\tiny},',...
        'legend style={font=\fontsize{7}{30}\selectfont, anchor=south, legend columns = 3, at={(0.4,1.03)}}']);
end