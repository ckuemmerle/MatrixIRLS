function plot_MC_v1(value_mat,parameters,alg_names,option,colors,problem,yname)

fdnames = fieldnames(parameters);
para1           =   fdnames{1};
values1         =   parameters.(para1);
nr_algos=length(alg_names);    

colorscheme= hsv(7);
if isempty(colors)
    colors=cell(1,20);
    colors{1}='k';      lines{1}='-';
    colors{2}='r';      lines{2}='-';
    colors{3}='k';      lines{3}='--';
    colors{4}='r';      lines{4}='--';
    colors{5}='k';      lines{5}=':';
    colors{6}='r';      lines{6}=':';
    colors{7}='r';      lines{7}=':';
    colors{8}='m';        lines{8}='-';
    colors{9}=colorscheme(2,:);         lines{9}='-';
    colors{10}=colorscheme(3,:);        lines{10}='-';
    colors{11}=[255/255,165/255,0/255];lines{11}='-'; 
    colors{12}=colorscheme(4,:);  lines{12}='-'; 
    colors{13}=colorscheme(6,:);        lines{13}='-';
    colors{14}=colorscheme(7,:);        lines{14}='-';
    colors{15}='o--m';
    colors{16}='d-.k';
    colors{17}='d-.r';
    colors{18}='d-.g';
    colors{19}='d-.b';
    colors{20}='d-.m';  
end
 figure('Name',['Success rates for and varying values of ',para1]);
    MaximizeFigureWindow();  
for i=1:nr_algos
if strcmp(option,'logy')
    semilogy(values1,value_mat(i,:),'LineStyle',lines{i},'LineWidth',5);
else
    plot(values1,value_mat(i,:),'LineStyle',lines{i},'LineWidth',5);
end
    hold on;
end
titlestr=['d1=',num2str(problem.d1),', d2=',num2str(problem.d1),...
    ', r=',num2str(problem.r),', kappa=',num2str(problem.cond_nr)];
title(titlestr)
[~,hObj]=legend(alg_names,'Interpreter','Latex');
 xlabel(para1,'interpreter','tex');
 ylabel(yname,'interpreter','tex');
 set(gcf,'name',titlestr) 

end     


