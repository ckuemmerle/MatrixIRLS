% This script plots the experimental data of Figure 1 contained in the 
% paper "Escaping Saddle Points in Ill-Conditioned Matrix Completion with 
% a Scalable Second Order Method" by Christian Kuemmerle and Claudio M. 
% Verdun, published at the ICML 2020 Workshop on 'Beyond First Order Methods in 
% Machine Learning'.
%
% In the experiment, the performance of several algorithms for exact matrix
% completion of (1000 x 1000) matrices of rank 5 with condition number
% kappa = 10 is compared in terms of median rel. Frobenius error of
% reconstruction vs. oversampling factor rho (which is the quotient between
% sample size and number of degrees of freedom).
%
% To re-run the experiments, see file
% "run_experiment_ICML2020_Figure1.m".
%
% Author: Christian Kuemmerle, June 2020.

clear all;
custom_options = struct;
custom_options.ColorOrderIndices=[2,3,4,5,6,7,8,6];

load(['experiment_MCalgos_ICML2020_Fig1_samplesvsmedian_kappa10.mat']);
plot_results_flexible(median_err_rel,parameters,alg_names,'logy',...
['Median of rel. Frob. errors of $\mathbf{X}^{(K)}$'],...
'Oversampling factor $\rho$',custom_options);
set(gcf,'Position',[744,820,384,230]);

save_filename=strcat('experiment_MCalgos_ICML2020_Fig1_samplesvsmedian_kappa10');
savefig(strcat(save_filename,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat(save_filename,'.tex'),...
        'height', '\figureheight', 'width', '\figurewidth',...
        'extraaxisoptions',['xlabel style={font=\tiny},',...
        'ylabel style={font=\tiny},',...
        'legend style={font=\fontsize{7}{30}\selectfont, anchor=south, legend columns = 4, at={(0.4,1.03)}}']);
end