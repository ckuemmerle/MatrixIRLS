function postprocess_fig(alg_names,xname,yname,varargin)

% We have two options:  option == '1d' means that we plot, in the case that
% we want to vary para2, too, several 1d plots.
%                       option == '2d' means that we plot one 2d plot in
%                       this case.


% load(parameterfile); % loads parameters

% % plot results
% instancesize    =   parameters.instancesize;
% d1               =   parameters.d1;
% d2           =   parameters.d2;
% m           =   parameters.m;
% r               =   parameters.r;
% tol             =   parameters.tol;
% N0              =   parameters.N0;
% 
% 
% 
markers = {'-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v', ...
           '-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v'};
 
if nargin >= 4 
  custom_options = varargin{1};
  if isfield(custom_options,'markers')
      markers =  custom_options.markers;
  end
  if isfield(custom_options,'ColorOrderIndices')
      ColorOrderIndices =  custom_options.ColorOrderIndices;
  end
else
    custom_options = [];
end


colorscheme = [0.00000 0.44700 0.74100
                   0.85000 0.32500 0.09800
                   0.92900 0.69400 0.12500
                   0.49400 0.18400 0.55600
                   0.46600 0.67400 0.18800
                   0.30100 0.74500 0.93300
                   0.63500 0.07800 0.18400%ScaledGD
                   0.08000 0.39200 0.25100
                   0.00000 0.44700 0.74100
                   0.85000 0.32500 0.09800
                   0.92900 0.69400 0.12500
                   0.49400 0.18400 0.55600];%
set(groot,'defaultAxesColorOrder',colorscheme)
% fig=gcf;
nr_algos = length(alg_names);
%         successrate_matrix = zeros(1,length(values1),nr_algos);
ax = gca;
for l=1:nr_algos
    if isfield(custom_options,'ColorOrderIndices')
        ax.Children(nr_algos+1-l).Color = colorscheme(l,:);
    end
    if isfield(custom_options,'markerssize')
        ax.Children(nr_algos+1-l).MarkerSize = custom_options.markerssize{l};
    end
end
%     titlestr=['d1=',num2str(problem.d1),', d2=',num2str(problem.d1),...
%         ', r=',num2str(problem.r),', kappa=',num2str(problem.cond_nr)];
%     title(titlestr)
%     [~,hObj]=legend(alg_names,'Interpreter','Latex');%,'FontSize',19);
xlabel(xname,'interpreter','Latex');
ylabel(yname,'interpreter','Latex');

end     
