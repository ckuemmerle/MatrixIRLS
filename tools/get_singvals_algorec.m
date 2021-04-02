function [sings,singsX0] = get_singvals_algorec(Xr,X0,sps_plc,r,iter,...
                                plotflag,varargin)
% This function calculates r singular values of matrix estimates given by
% the iteration number 'iter' of cell array Xr, and of the matrix X0. If
% plotflag is true, provide a scatter plot of the singular values.
% =========================================================================
% Parameters
% ----------
% Returns
% ----------
% =========================================================================
% Author: Christian Kuemmerle, 2020.
nr_algos = length(Xr);

if isfloat(iter)
    for l=1:nr_algos
        if length(Xr{l}) < iter
            error('Choose a smaller iteration number.')
        end
    end
end
sings=cell(1,nr_algos);
for l=1:nr_algos
    if isfloat(iter)
        sings{l} = get_singvals(Xr{l}{iter},sps_plc,r);
    else
        if strcmp(iter,'last')
            sings{l} = get_singvals(Xr{l}{end},sps_plc,r);
        elseif strcmp(iter,'first')
            sings{l} = get_singvals(Xr{l}{1},sps_plc,r);
        end
    end
end
singsX0 = get_singvals(X0,sps_plc,r);

if plotflag
    figure;hold on
    for l=1:nr_algos
        s = scatter(1:r,log10(sings{l}(1:r)),40);%,'LineWidth',5);
%         s.LineWidth = 5;
        hold on;
    end
    scatter(1:r,log10(singsX0(1:r)),100,'x');%,'x','LineWidth',10
    xlabel('Index of singular value');
    ylabel('Log. of value of singular value');
    if nargin >= 7
        alg_names = varargin{1};
        legend(alg_names,'FontSize',12,'Interpreter','TeX');
    end
%     legend(alg_name,'FontSize',12,'Interpreter','TeX'); %,'Interpreter','LaTeX'
    title('First singular values of reconstuctions')
    set(gcf,'name','First singular values of reconstuctions') 
        
    figure;hold on
    for l=1:nr_algos
        scatter(1:r,log10(abs(sings{l}(1:r)-singsX0(1:r))./singsX0(1:r)),40);%,'LineWidth',5);
%         s.LineWidth = 5;
        hold on;
    end
    xlabel('Index of singular value');
    ylabel('Rel. error to singular value of X0');
    if nargin >= 7
        alg_names = varargin{1};
        legend(alg_names,'FontSize',12,'Interpreter','TeX');
    end
%     legend(alg_name,'FontSize',12,'Interpreter','TeX'); %,'Interpreter','LaTeX'
    title('Rel. errors to singular value of X0')
    set(gcf,'name','Rel. errors to singular value of X0') 
end

end

