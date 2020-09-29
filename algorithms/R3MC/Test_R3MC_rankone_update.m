% This is a generic test file for calling algorithm R3MC
% for the low-rank matrix completion problem.
%
% Refer "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
% B. Mishra and R. Sepulchre,
% Technical report, arXiv:1306.2672, 2013.
%
% Author:
% Bamdev Mishra
% b.mishra@ulg.ac.be



clear all; close all; clc;
rndstr = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(rndstr);  % fixed seed




%% Problem specifications

d1 = 5000;
d2 = 5000;
true_rank = 10;
over_sampling = 4;
CN = 10^7 ; 200; 300; 500; % Condition number

fprintf('Rank %i matrix of size %i times %i and OS = %i\n', true_rank, d1, d2, over_sampling);

tol = 1e-18; % Tolerance; Look at the R3MC file to stop R3MC in a differente way
verbosity = true; % Show output
max_iter = 1000; % Maximum number of iterations allowed
ls_maxiter = 100; % Maximum linesearch allowed at each iteration
beta_type = 'P-R'; % Algorithm type. Other choices include H-S, F-R and off
linearized_linesearch = true; % Whether to use linearized stepsize search
accel_linesearch = true; % Whether to use accelerated linesearch
compute_predictions = true; % If interested in compute recovery at each iteration
random_initialization = false; % How to initialize the iterates
ORTH_VALUE = Inf; % Inf gives the best results for R3MC
sigma_armijo = 1e-4; % for Armijo linesearch. Keep a low valaue while going for linearized stpesize search
if ~linearized_linesearch,
    sigma_armijo = 0.5;
end




%% Generate data well-conditioned or ill-conditioned

dof = (d1 + d2 - true_rank)*true_rank;
prop_known = over_sampling * dof / (d1 * d2);

% Creating a well/ill-conditioned matrix
prob_options.type ='decaying'; % Ill-coniditioned matrix depending on CN

if strcmp(prob_options.type,'decaying'),    
    S0 = 1000*diag(logspace(-log10(CN),0,true_rank)); fprintf('Exponential decay of singular values.\n'); % Exponential decay
%         S0 = diag(1:(CN-1)/(true_rank -1):CN); fprintf('Linear decay of singular values.\n'); % Linear decay
    
    prob_options.diag = S0;
    prob_options.scaling = 1;
    S0 = prob_options.scaling*prob_options.diag;
    fprintf('Creating a matrix with singular values...\n')
    for kk = 1: length(diag(S0));
        fprintf('%s \n', num2str(S0(kk, kk), '%10.5e') );
    end
end

m = round(prop_known*d1*d2);
Matcom = prob_matrix(d1,d2,true_rank,m,prob_options); % Problem related structure variable
singular_vals = svd(Matcom.left'*Matcom.left);
condition_number = sqrt(max(singular_vals)/min(singular_vals));
fprintf('Condition number is %f \n', condition_number);




%% Sparse structure

sparse_structure = sparse(Matcom.row', Matcom.col', Matcom.values_Omega, d1, d2, m);





%% Initialization





%% Interface with our algorithms

data_ls.rows = Matcom.row';
data_ls.cols = Matcom.col';
data_ls.entries = Matcom.values_Omega';
data_ls.nentries = length(data_ls.entries);

% Testing data
data_ts.nentries = 1*data_ls.nentries;
data_ts.rows = randi(d1, data_ts.nentries, 1);
data_ts.cols = randi(d2, data_ts.nentries, 1);
data_ts.entries = partXY(Matcom.left', Matcom.right',data_ts.rows,data_ts.cols,data_ts.nentries)';





%% R3MC

% Mandatory parameters
model.d1 = d1; % Number of rows
model.d2 = d2; % Number of columns
% model.r = r; % Rank 
% model.U = U;
% model.V = V;
% model.R = B;
model.M = sparse_structure; % Sparse structure


% Optional parameters
params.tol = tol;
params.maxiter = max_iter;
params.sigma_armijo = sigma_armijo;
params.compute_predictions = compute_predictions;
params.beta_type = beta_type;
params.linearized_linesearch =  linearized_linesearch;
params.verb = verbosity;
params.orth_value = ORTH_VALUE;
params.accel_linesearch = true;
params.ls_maxiter = ls_maxiter;
params.gamma = 0;

% % Call R3MC
% fprintf('--------------- R3MC -------------\n');
% tic;
% [model_R3MC, infos_R3MC] = R3MC(data_ts, data_ls, model, params);
% toc;

rank_max =true_rank;


test_error_rank_urv = [];
ranklist_urv = [];

for jj = 1 : rank_max
    fprintf('>>> Rank %i\n', jj);
    
    if jj == 1
        r = jj;

        if ~random_initialization
            [U,B,V] = svds(sparse_structure, r);
            G = U*(B.^(0.5));
            H = V*(B.^(0.5));
            fprintf('**** Initialization by taking %i dominant SVD\n', r);
        else
            G = randn(d1, r); H = randn(d2, r);
            [Qg,Rg] = qr(G, 0);
            [Qh,Rh] = qr(H, 0);
            [q1,b,q2] = svd(Rg * Rh');
            U = Qg * q1;
            V = Qh * q2;
            B = b;
            fprintf('**** Random initialization\n');
        end
        model.U = U;
        model.V = V;
        model.R = B;
    end
    model.r =jj;
    tic;
    [model_R3MC, infos_R3MC] = R3MC(data_ts, data_ls, model, params);
    toc;
    
    
    % Collect statistics at every rank
%     test_error_rank_urv = [test_error_rank_urv; infos_CG_urv.test_error(end)];
    ranklist_urv = [ranklist_urv ; size(model.R, 1)];

    % Print
    fprintf('Test error %.5e \n', infos_R3MC.test_error(end));


    model = model_R3MC;
    
    % Rank update
    [u, sig, v] = svds(model_R3MC.M, 1);

    [model.U, model.R, model.V] = svd_rank_1_update(model.U, model.R, model.V, -sig*u, v); % Rank update along the negative gradient direction

    clear u sig v
    
end


%% Plots

% Cost versus time
fs = 20;
figure;
semilogy((infos_R3MC.iter_time),(infos_R3MC.costs),'-','Color','blue','LineWidth',2);
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel(ax1,'Time in seconds','FontSize',fs);
ylabel(ax1,'Cost','FontSize',fs);
axis([get(gca,'XLim') tol 1e3])
legend('R3MC');
legend 'boxoff';
title([num2str(d1),' by ',num2str(d2),', rank ',num2str(true_rank),', OS ',num2str(over_sampling),', Condition No. ', num2str(condition_number)])
box off;


% Cost versus iterations
fs = 20;
figure;
semilogy(0:length(infos_R3MC.costs)-1,(infos_R3MC.costs),'-','Color','blue','LineWidth',2);
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel(ax1,'Number of iterations','FontSize',fs);
ylabel(ax1,'Cost ','FontSize',fs);
axis([get(gca,'XLim') tol 1e3])
legend('R3MC');
legend 'boxoff';
box off;
title([num2str(d1),' by ',num2str(d2),', rank ',num2str(true_rank),', OS ',num2str(over_sampling),', Condition No. ', num2str(condition_number)])


% Test error versus iterations
if compute_predictions
    fs = 20;
    figure;
    semilogy(0:length(infos_R3MC.test_error)-1,(infos_R3MC.test_error),'-','Color','blue','LineWidth',2);
    ax1 = gca;
    set(ax1,'FontSize',fs);
    xlabel(ax1,'Number of iterations','FontSize',fs);
    ylabel(ax1,'Test error','FontSize',fs);
    legend('R3MC');
    legend 'boxoff';
    box off;
    title([num2str(d1),' by ',num2str(d2),', rank ',num2str(true_rank),', OS ',num2str(over_sampling),', Condition No. ', num2str(condition_number)])
    
    
    
    fs = 20;
    figure;
    semilogy((infos_R3MC.iter_time),(infos_R3MC.test_error),'-','Color','blue','LineWidth',2);
    ax1 = gca;
    set(ax1,'FontSize',fs);
    xlabel(ax1,'Time in seconds','FontSize',fs);
    ylabel(ax1,'Test error','FontSize',fs);
    legend('R3MC');
    legend 'boxoff';
    title([num2str(d1),' by ',num2str(d2),', rank ',num2str(true_rank),', OS ',num2str(over_sampling),', Condition No. ', num2str(condition_number)])
    box off;
end




