% This script reproduces the experiment of Figure 4 of the paper
% C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples", ICML 2021.
%
% It performs a comparison of the execution time of different algorithms 
% (in seconds) as a function  of the matrix size. In the paper, we compare
% R2RILS and MatrixIRLS for matrices of size m x (m+100)  with rank=5 or 
% rank=10 and condition number kappa=1e2. Every point represents the
% average of 50 experiments.

% It first runs MatrixIRLS and saves the corresponding information. Then it
% does the same with R2RILS and, in the end, it plots the two algorithms 
% together for rank=5 and rank=10.
%
% Author: Claudio M. Verdun, June 2020.
%% Choose parameters for the algorithm 'MatrixIRLS'
alg_names={'MatrixIRLS'};
opts_custom.p=0; %%% choose (non-)convexity parameters p for IRLS
p=opts_custom.p;
opts_custom.tol = 1e-6;
opts_custom.N0 = 500;
opts_custom.N0_firstorder = 4000;
opts_custom.N0_inner = 150; 
opts_custom.mode_linsolve = 'tangspace';
opts_custom.type_mean={'geometric'};
opts_custom.epsmin = 1e-16;
opts_custom.lambda = 0;
opts_custom.saveiterates = 0;
opts_custom.verbose = 2;
%% Run MatrixIRLS for matrix completion
no_rows_cell={};
evaluation_time={};
evaluation_time_average={};
conv_experiments=[];

modeX0      = 'condition_control_linear';
complexflag = 0;
cond_nr = 1e2;
max_nr_resample = 1000;
verbose_errors = 1;
frob_mode = 'full';
rng('shuffle');
for rank=1:2
    for d1=100:100:1000
        for no_experiments=1:50
            d2=d1+100;
            no_rows_cell={no_rows_cell{:},d1};
            r = 5*rank;  % Rank:
            k = 2;
            df_LR = @(rr) rr*(d1 + d2 - rr);  % Number of degrees of freedom of the setting:
            df_LR_val=df_LR(r);
            m = floor(min(2.5*df_LR_val,d1*d2)); % Number of measurements:
            [Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
            [rowind,colind] = ind2sub([d1,d2],Omega);
            % sample the ground truth matrix X0 and measured entries y
            opts_custom.tol_CG_fac=1e-5*cond_nr^(-1);
            [U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
            X0 = {U0,V0};
            y = partXY(U0',V0',rowind,colind,m).';
            [Xr,outs] = run_MC_algos(Phi,y,r,alg_names,opts_custom);
            % sum the time of all iterations and concatenate them into a cell: 
            evaluation_time={evaluation_time{:},outs{1,1}.time(end)};
            [error_fro_rel,error_fro] = get_frob_errors(Xr,X0,Phi,alg_names,...
                frob_mode,verbose_errors);
            counter_d1=d1/100;
            conv_experiments(rank,counter_d1,no_experiments)=error_fro_rel{1}(end);
        end
        time_matrix=cell2mat(evaluation_time);
        time_average_matrix=mean(time_matrix);
        evaluation_time_average={evaluation_time_average{:},time_average_matrix};
    end
end   


conv_experiments_MatrixIRLS = conv_experiments;
evaluation_time_average_MatrixIRLS=evaluation_time_average;
time_MatrixIRLS=cell2mat(evaluation_time_average_MatrixIRLS);
time_MatrixIRLS_rank5=time_MatrixIRLS(1:10);
time_MatrixIRLS_rank10=time_MatrixIRLS(11:20);

%% Run R2RILS for matrix completion

alg_names={'R2RILS'};
no_rows_cell={};
evaluation_time={};
evaluation_time_average={};
conv_experiments=[];
for rank=1:2
    for d1=100:100:1000
        for no_experiments=1:50
            d2=d1+100;
            no_rows_cell={no_rows_cell{:},d1};
            r = 5*rank;  % Rank:
            df_LR = @(rr) rr*(d1 + d2 - rr);  % Number of degrees of freedom of the setting:
            df_LR_val=df_LR(r);
            m = floor(min(2.5*df_LR_val,d1*d2)); % Number of measurements:
            [Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
            [rowind,colind] = ind2sub([d1,d2],Omega);
            % Sample the ground truth matrix X0 and measured entries y
            opts_custom.tol_CG_fac=1e-5*cond_nr^(-1);
            [U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
            X0 = {U0,V0};
            y = partXY(U0',V0',rowind,colind,m).';
            [rowind,colind] = find(Phi);
            [Xr,outs] = run_MC_algos(Phi,y,r,alg_names,opts_custom);
            evaluation_time={evaluation_time{:},outs{1,1}.time(end)};
            [error_fro_rel,error_fro] =get_frob_errors(Xr,X0,Phi,alg_names,...
                frob_mode,verbose_errors);
            counter_d1=d1/100;
            conv_experiments(rank,counter_d1,no_experiments)=error_fro_rel{1}(end);
        end
        time_matrix=cell2mat(evaluation_time);
        time_average_matrix=mean(time_matrix);
        evaluation_time_average={evaluation_time_average{:},time_average_matrix};
    end
end   

conv_experiments_R2RILS = conv_experiments;
evaluation_time_average_R2RILS=evaluation_time_average;
time_R2RILS=cell2mat(evaluation_time_average_R2RILS);
time_R2RILS_rank5=time_R2RILS(1:10);
time_R2RILS_rank10=time_R2RILS(11:20);

x=100:100:1000;
figure
plot(x,time_MatrixIRLS_rank5,'-o', 'LineWidth',1.5,'Color',[0 0.44 0.74])
hold on
plot(x,time_R2RILS_rank5,'-o','LineWidth',1.5,'Color',[0.85 0.32 0.098])
hold on
plot(x,time_MatrixIRLS_rank10,'LineWidth',2,'Color',[0 0.44 0.74])
hold on
plot(x,time_R2RILS_rank10,'LineWidth',2,'Color',[0.85 0.32 0.098])
grid on
legend('MatrixIRLS r=5','R2RILS r=5','MatrixIRLS r=10','R2RILS r=10')
title('Execution Time of Different Algorithms for \kappa=100')
xlabel('m (with n=m+100)') 
ylabel('Execution time (seconds)')

%% Save data and figure
outs = remove_intermediate_iterates(outs);
clear Xr;
filename=strcat('experiment_MCalgos_ICML2021_Fig4_compareMIRLS_R2RILS_rk5and10');
save(strcat('results/',filename,'.mat'))
savefig(strcat('results/',filename,'.fig'))
if exist('matlab2tikz')
    matlab2tikz(strcat('results/',filename,'.tex'))
end