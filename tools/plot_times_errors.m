function fig = plot_times_errors(outs,error_l2_rel,alg_name,varargin)
%plot_times_errors

%% Get options
if nargin >= 4
    plot_options = varargin{1};
else
    plot_options = [];
end
if isfield(plot_options,'markers')
    markers = plot_options.markers;
else
    markers = {'-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v', ...
       '-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v'};
end
if isfield(plot_options,'markersize')
    markersize = plot_options.markersize;
else
    markersize = 12;
end
if isfield(plot_options,'maxtime')
    maxtime = plot_options.maxtime;
else
    maxtime = false;
end

%% Determine which algorithms have time stamps
nr_algos=length(alg_name);
times=cell(1,nr_algos);
plot_flag = zeros(1,nr_algos);
for l=1:nr_algos
    if isfield(outs{l},'time') && (length(outs{l}.time) == length(error_l2_rel{l}))
        times{l}=outs{l}.time;
        plot_flag(l) = 1;
    elseif isfield(outs{l},'time') && (length(error_l2_rel{l}) == 1)
        times{l}=outs{l}.time(end);
        plot_flag(l) = 1;
    else
        times{l}=[];
        plot_flag(l) = 0;
    end
end
%% Plot errors of algorithms with time stamps
%colors = {[0 0 1], [0 0 0.5],[0 0 0],[1 0 0],[0 1 0]};
fig = figure;
alg_name_plot = {};
subplot('Position',[0.15 0.2 0.8 0.7]);
for l = 1:nr_algos
    if plot_flag(l)
        semilogy(times{l},error_l2_rel{l},markers{l},'MarkerSize',markersize);hold on;
        alg_name_plot = {alg_name_plot{:},alg_name{l}};
    end
end
xlabel('Time in seconds');
ylabel('Relative $\ell_2$-error','Interpreter','latex');
legend(alg_name_plot,'Interpreter','latex','FontSize',12);
title('Relative error vs. time')
set(gcf,'name','Rel. error vs. time')
if maxtime
    xlim([0 maxtime])
end
end

