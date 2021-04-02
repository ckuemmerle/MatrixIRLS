function visualize_errorcurves_combined(error_fro_rel,alg_name,varargin)

if nargin == 3
   problem = 'RobustPCA';
   nr_outliers = varargin{1};
elseif nargin > 3
   problem = [];
   titlestring = varargin{2};
else
   problem = [];
   titlestring = 'Relative error vs. iteration count';
end

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
nr_algos=length(alg_name);
if ~isequal(nr_algos,length(error_fro_rel))
   error('Mismatch between provided number of algorithms and number of calculated error statistics.');
end

N=cell(1,nr_algos);
for i = 1:nr_algos
    N{i}=length(error_fro_rel{i});
end

maxIt=N{1};
minIt=N{1};
for l=1:nr_algos
   maxIt=max(maxIt,N{l});
   minIt=min(minIt,N{l});
end
markers = {'-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v', ...
           '-x', '-+', '-*', '-o','-x', '-s', '-d', '-^', '-v'};

% bigFig;
figure;

if strcmp(problem,'RobustPCA')
    subplot(1,2,1)
    for l=1:nr_algos
       semilogy(error_fro_rel{l},markers{l},'MarkerSize',8);
       hold on
    end
    xlabel('iterations');
    ylabel('Logarithm of relative Frobenius error');
    set(gca,'xtick',0:floor(maxIt./15):maxIt,'xticklabel',0:floor(maxIt./15):maxIt);
    legend(alg_name,'FontSize',12,'Interpreter','TeX');
    title('Relative error vs. iterations')
    subplot(1,2,2)
    for l=1:nr_algos
       plot(nr_outliers{l},markers{l},'MarkerSize',8);
       hold on
    end
    xlabel('iterations');
    ylabel('Number of outliers');
    set(gca,'xtick',0:floor(maxIt./15):maxIt,'xticklabel',0:floor(maxIt./15):maxIt);
    legend(alg_name,'FontSize',12);
    title('Nr. of outliers vs. iterations')
        set(gcf,'name','Algorithmic comparisons for Robust PCA') 
else
%     subplot(1,2,1)
%     for l=1:nr_algos
%         plot(error_fro_rel{l},markers{l},'MarkerSize',8);
%         hold on
%     end
%     xlabel('iterations');
%     ylabel('Relative Frobenius error');
%     set(gca,'xtick',0:floor(maxIt./15):maxIt,'xticklabel',0:floor(maxIt./15):maxIt);
%     legend(alg_name,'FontSize',12); %,'Interpreter','LaTeX'
%     title('Relative error vs. iteration count')

%     subplot(1,2,2)
    for l=1:nr_algos
       semilogy(error_fro_rel{l},markers{l},'MarkerSize',8); %
       hold on
    end
    xlabel('iterations');
    ylabel('Relative Frobenius error (log-scale)');
    set(gca,'xtick',0:floor(maxIt./15):maxIt,'xticklabel',0:floor(maxIt./15):maxIt);
    legend(alg_name,'FontSize',12,'Interpreter','TeX'); %,'Interpreter','LaTeX'
    title(titlestring)
    set(gcf,'name',titlestring) 
end
        
end

% function plot_logaplot
% 
% end