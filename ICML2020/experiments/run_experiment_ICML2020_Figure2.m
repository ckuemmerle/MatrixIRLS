% This script repeats the experiment illustrated in Figure 1 of the paper
% "Escaping Saddle Points in Ill-Conditioned Matrix Completion with 
% a Scalable Second Order Method" by Christian Kuemmerle and Claudio M. 
% Verdun, published at the ICML 2020 Workshop on 'Beyond First Order Methods in 
% Machine Learning'.
%
% In the experiment, the performance of several algorithms for exact matrix
% completion of (1000 x 1000) matrices of rank 5 with condition number
% kappa = 10^5 is compared in terms of median rel. Frobenius error of
% reconstruction vs. oversampling factor rho (which is the quotient between
% sample size and number of degrees of freedom).
%
% Caution: The execution of this experiment may take many hours! It uses
% the parallel computing toolbox, and runtime will depend on the number of
% available cores.
% See also: File experiment_ICML2020_Figure1_samplesvsmedian_kappa100000.m
%
% Author: Christian Kuemmerle, June 2020.
clear all;
%%% Define how many jobs run in parallel
pc=parcluster('local');
n_jobs=pc.NumWorkers; 
delete(gcp('nocreate'));

%%% Choose problem parameters and run experiment
instancesize = 50; 
alg_names = {'MatrixIRLS','R2RILS','RTRMC','LMaFit','ScaledASD','ScaledGD','NIHT'}; 
problem =struct; 
problem.d1 = 1000; 
problem.d2 = 1000; 
problem.r = 5; 
problem.modeX0 = 'condition_control_log'; 
problem.cond_nr = 10^5; 
parameters.rho=[1.05,1.1:0.1:4.0]; 
filename = MC_samplecomp_experiment(n_jobs,instancesize,alg_names,parameters,problem);


load(fullfile('results',strcat(filename,'.mat')));
save('experiment_MCalgos_ICML2020_Fig2_samplesvsmedian_kappa100000.mat')