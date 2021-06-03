function plot_results_flexible(value_mat,parameters,alg_names,option,yname,varargin)

markers = {'-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v', ...
           '-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v'};

if nargin >= 6
   xname = varargin{1}; 
   if nargin >= 7 
      custom_options = varargin{2};
      if isfield(custom_options,'markers')
          markers =  custom_options.markers;
      end
      if isfield(custom_options,'ColorOrderIndices')
          ColorOrderIndices =  custom_options.ColorOrderIndices;
      end
   end
end


fdnames = fieldnames(parameters);
para1           =   fdnames{1};
values1         =   parameters.(para1);
nr_algos=length(alg_names);     
colorscheme = [0.00000 0.44700 0.74100
               0.85000 0.32500 0.09800
               0.92900 0.69400 0.12500
               0.49400 0.18400 0.55600
               0.46600 0.67400 0.18800
               0.30100 0.74500 0.93300
               0.63500 0.07800 0.18400
               0.08000 0.39200 0.25100
               0.00000 0.00000 0.00000];
set(groot,'defaultAxesColorOrder',colorscheme)
figure 
for i=1:nr_algos
    if strcmp(option,'logy')
        semilogy(values1,value_mat(i,:),markers{i},'LineWidth',1);
    else
        plot(values1,value_mat(i,:),markers{i},'LineWidth',1);
    end
    hold on;
    if isfield(custom_options,'ColorOrderIndices')
        set(gca,'ColorOrderIndex',ColorOrderIndices(i));
    end
end

[~,hObj]=legend(alg_names,'Interpreter','Latex');
if nargin >= 6
    xlabel(xname,'interpreter','Latex');
else
     xlabel(para1,'interpreter','Latex');
end
     ylabel(yname,'interpreter','Latex');

end     
