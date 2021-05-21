function [err_fro_rel,erPhi_fro_rel,erPhic_fro_rel,...
    err_tail_rel,times]...
    = MC_run_pipeline(alg_names,alg_opts,...
    problem,instancesize,n_jobs,varargin)
% Runs the algorithms indicated by the cell array 'alg_names'
% for the problem type described by the struct 'problem'. Of this problem
% type, a number of 'instancesize' independent instances are executed.
% Changes from the standard algorithmic parameters are contained in
% 'alg_opts'. The max. number of parallel jobs to be used is 'n_jobs'.
% =========================================================================
% Parameters
% ----------
% alg_names:    cell array, size (1 x nr_algos). Each cell contains a 
%               character string with the name of an algorithm to be used.
% alg_opts:     struct. Fields indicate algorithmic options to override
%               the standard ones (for any algorithm).
% problem:  struct. Contains problem parameters, fields include:
%           - d1    int. First dimension of matrix to completed
%           - d2    int. Second dimension of matrix to completed
%           - r     int. Rank of matrix to be completed
%           - modeX0 char string or float/int. Indicates the random model
%                   for matrix to be completed.
%           - cond_nr. float. Fixes the condition number of matrix to be
%                   completed (only used if 'modeX0' is a character string).
% instancesize: int. Number of independent random instance to be run in the
%               experiment.
% n_jobs:       int. Max. number of parallel jobs to be used in the
%               parallel execution of the experiment (parallelization is
%               used across the different instances).
% Optional argument: struct with a field that contains an array. The
%               different values of this field indicate a problem parameter
%               that is varied in the experiment.
% Returns
% ----------
% err_fro_rel:  cell array, size (nr_algos x nr_parameters x instancesize).
%               'nr_parameters' is the length of the field of optinal
%               argument. 
%               Contains the relative Frobenius error to ground truth 
%               matrix of reconstruction in its (i,j,k)-th entry 
%               corresponding to the i-th algorithm, j-th parameter and 
%               k-th instance.
% erPhi_fro_rel: cell array. As 'err_fro_rel', but corresponds to relative
%               Frobenius errors restricted to the known indices 'Omega'.
% erPhic_fro_rel: cell array. As 'err_fro_rel', but corresponds to relative
%               Frobenius errors restricted to the unknown indices 'Omega^c'.
% err_tail_rel: cell array. As 'err_fro_rel', but corresponds to relative 
%               'tail Frobenius errors' of order floor(r/2), cf. function
%               'get_tail_frob_errors'.
% times:        cell array, size (nr_algos x nr_parameters x instancesize).
%               Reports the times until convergence of algorithm at
%               respective problem instance.
% =========================================================================
% Author: Christian Kuemmerle, 2020.

nr_algos = length(alg_names);

if nargin >= 6
    if ~isstruct(varargin{1})
        error('Input variable must be a struct (containing the options to be changed).')
    end
    parameters=varargin{1};
    para_names=fieldnames(parameters);
    if length(para_names) > 1
        error('Please specify only at most one hyperparameter.')
    end
end

err_fro_rel = cell(nr_algos,length(parameters.(para_names{1})),...
                    instancesize);
erPhi_fro_rel = cell(nr_algos,length(parameters.(para_names{1})),...
                    instancesize);
erPhic_fro_rel = cell(nr_algos,length(parameters.(para_names{1})),...
    instancesize);

err_tail_rel = [];
% err_tail_rel = cell(nr_algos,length(parameters.(para_names{1})),...
%     instancesize);
times = cell(nr_algos,length(parameters.(para_names{1})),...
    instancesize);

delete(gcp('nocreate'));
pc = parcluster('local');
% set number of workers
pc.NumWorkers = n_jobs; %nw;
% start the parallel pool
poolobj = parpool(pc, n_jobs); %nw

for j=1:length(parameters.(para_names{1}))
    if length(para_names) == 1
        problem = overwrite_options_LowRank(problem,parameters,j,[]);
        disp(['Run ',num2str(instancesize),' instances of matrix completion for ',para_names{1},'=',...
            num2str(parameters.(para_names{1})(j)),'...'])
        r = problem.r;
        k_tail = max(1,ceil(r/2));
        d1 = problem.d1; 
        d2 = problem.d2;
        if strcmp(para_names{1},'rho')
            df_LR = @(rr) rr*(problem.d1 + problem.d2 - rr);
            m = round(problem.rho*df_LR(problem.r));
        else
            m = problem.m;
        end
        max_nr_resample = 1000;
        modeX0 = problem.modeX0;
        complexflag = 0;
        if strcmp(modeX0,'condition_control_1/x2') || strcmp(modeX0,'condition_control_linear') ...
            || strcmp(modeX0,'condition_control_log')
            cond_nr = problem.cond_nr;
            alg_opts.tol_CG_fac = alg_opts.tol_CG_fac./cond_nr;
        else
            cond_nr  = [];
            alg_opts.tol_CG_fac = alg_opts.tol_CG_fac;
        end
        parfor i=1:instancesize
            if strcmp(modeX0,'condition_control_1/x2') || strcmp(modeX0,'condition_control_linear') ...
                || strcmp(modeX0,'condition_control_log')
                [U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag,cond_nr);
            else
                [U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag);
            end
            X0 = {U0,V0};
            [Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
            [rowind,colind] = ind2sub([d1,d2],Omega);
            y = partXY(U0',V0',rowind,colind,m).';

            [Xr_c,outs_c] = run_MC_algos(Phi,y,r,alg_names,alg_opts);
            
            verbose = 0; % provide text output after error calculations
            error_fro_rel_c =get_frob_errors(Xr_c,X0,Phi,alg_names,'full',verbose);
            erPhi_fro_rel_c =get_frob_errors(Xr_c,X0,Phi,alg_names,'Phi',verbose);
            erPhic_fro_rel_c =get_frob_errors(Xr_c,X0,Phi,alg_names,'Phi_comp',verbose);
%             error_tail_rel_c =get_tail_frob_errors(Xr_c,X0,Phi,alg_names,...
%                 k_tail,'full',verbose);
            times_c = cell(nr_algos,1);
            for k=1:nr_algos
                if isempty(outs_c{k}.time)
                    times_c{k} = 0;
                else
                    times_c{k}=outs_c{k}.time(end);
                end
            end
            err_fro_rel(:,j,i) = error_fro_rel_c(:);
            erPhi_fro_rel(:,j,i) = erPhi_fro_rel_c(:);
            erPhic_fro_rel(:,j,i) = erPhic_fro_rel_c(:);
%             err_tail_rel(:,j,i) = error_tail_rel_c(:);
            times(:,j,i) = times_c(:);
        end
    end
end

% %===============================================================================
% % Close parallel pool
% %===============================================================================
delete(poolobj);
    
end

