function filename = MC_samplecomp_experiment(n_jobs,instancesize,...
    alg_names,parameters,problem)

alg_opts = getDefaultOpts_IRLS;
alg_opts.p = 0; %%% choose (non-)convexity parameters p for IRLS
alg_opts.tol=1e-9;
alg_opts.tol_CG_fac=1e-5;
alg_opts.N0= 400;
alg_opts.N0_inner=500;
alg_opts.N0_firstorder = 4000;
alg_opts.type_mean={'geometric'};
alg_opts.qmean_para= 1;%p/(p-2);
alg_opts.mode =  'objective_thesis';
alg_opts.mode_eps = 'oracle_model_order';
alg_opts.adaptive_cg = 0;
alg_opts.mode_linsolve = 'tangspace';
alg_opts.tracking = 0;
alg_opts.saveiterates = 0;
alg_opts.recsys = 0;
alg_opts.verbose = 0;

[err_fro_rel,erPhi_fro_rel,erPhic_fro_rel,...
    err_tail_rel,times]...
    = MC_run_pipeline(alg_names,alg_opts,...
    problem,instancesize,n_jobs,parameters);

%% Save data
curdate = datestr(now,'yyyy-mm-dd_HH-MM-SS'); % get string with current date and time
if isfield(problem,'cond_nr') && not(isempty(problem.cond_nr))
    filename = strcat('MCsamplecomp_d1_',num2str(problem.d1),...
        '_d2_',num2str(problem.d2),'_r_',num2str(problem.r),'_kappa_',num2str(problem.cond_nr),'_',curdate);
else
    filename = strcat('MCsamplecomp_d1_',num2str(problem.d1),...
        '_d2_',num2str(problem.d2),'_r_',num2str(problem.r),'_',curdate);
end
para_names=fieldnames(parameters);
nr_para_names=length(para_names);
filename_extension = [];
%%% Read all the parameter fieldnames and values, add to file name
for i=1:nr_para_names
    filename_extension=strcat(filename_extension,'_',para_names{i});
    para_value_first = parameters.(para_names{i})(1);
    para_value_last =  parameters.(para_names{i})(end);
    filename_extension=strcat(filename_extension,num2str(para_value_first),...
       '-',num2str(para_value_last));
end

%% Postprocessing
success_threshold = 1e-3;
successrate = get_successrate(err_fro_rel,success_threshold);
successrate_tail = get_successrate(err_tail_rel,success_threshold);

av_err_fro_rel     = average_errors(err_fro_rel);
av_erPhi_fro_rel   = average_errors(erPhi_fro_rel);
av_erPhic_fro_rel  = average_errors(erPhic_fro_rel);
av_err_tail_rel    = average_errors(err_tail_rel);
failurerate=1-successrate;
failurerate_tail=1-successrate_tail;
median_err_rel = get_median_errors(err_fro_rel);
median_err_rel_tail = get_median_errors(err_tail_rel);

filename = strcat(filename,filename_extension);
save(fullfile('results',strcat(filename,'.mat')))
disp('Save file with results!')
plot_MC_v1(failurerate,parameters,alg_names,'1d',[],problem,...
    ['Failure rate w/ threshold ',num2str(success_threshold)]);
savefig(fullfile('results',strcat(filename,'_failurerate.fig')));
print(fullfile('results',strcat(filename,'_failurerate.eps')),'-depsc');
plot_MC_v1(failurerate_tail,parameters,alg_names,'1d',[],problem,...
    ['Failure rate (tail) w/ threshold ',num2str(success_threshold)]);
savefig(fullfile('results',strcat(filename,'_failurerate_tail.fig')));
print(fullfile('results',strcat(filename,'_failurerate_tail.eps')),'-depsc')
plot_MC_v1(median_err_rel,parameters,alg_names,'logy',[],problem,...
    ['Median of relative Frobenius errors']);
savefig(fullfile('results',strcat(filename,'_medianerr.fig')));
print(fullfile('results',strcat(filename,'_medianerr.eps')),'-depsc');

plot_MC_v1(median_err_rel_tail,parameters,alg_names,'logy',[],problem,...
    ['Median of relative tail Frobenius errors']);
savefig(fullfile('results',strcat(filename,'_medianerrtail.fig')));
print(fullfile('results',strcat(filename,'_medianerrtail.eps')),'-depsc');
disp('Figures saved successfully!')
end

function error_avg = average_errors(error_cell)
%average_errors Summary of this function goes here
%   Detailed explanation goes here
[nr_algos,parameter_nrs,instancesize]=size(error_cell);
error_avg=cell(nr_algos,parameter_nrs);
for k=1:nr_algos
    for j=1:parameter_nrs
        error_avg{k,j}=0;
        for i=1:instancesize
            error_avg{k,j}=error_avg{k,j}+error_cell{k,j,i};
        end
        error_avg{k,j}=error_avg{k,j}./instancesize;
    end
end
end

function successrate = get_successrate(error_cell,success_threshold)
[nr_algos,nr_paras,instancesize] = size(error_cell);
successrate = zeros(nr_algos,nr_paras);
for j=1:nr_paras
    succ = zeros(nr_algos,1);
    for k=1:nr_algos
        for i=1:instancesize
            if error_cell{k,j,i}(end) < success_threshold
                succ(k)=succ(k)+1;
            end 
        end
    end
    successrate(:,j)=succ./instancesize;
end
end