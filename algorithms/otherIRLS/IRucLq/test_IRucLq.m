%   copyright (c) by Yangyang Xu, Oct. 2012
%   yangyang.xu@rice.edu
%
% Paper:
% M.-J. Lai, Y. Xu, W. Yin. 
% Improved iteratively reweighted least squares for unconstrained smoothed lq
% minimization, SIAM Journal on Numerical Analysis, 5(2), 927-957, 2013.

clear; close all;

%% sparse vector recovery
m = 64; % sample size
n = 128; % vector size
k_seq = 20:2:32; % number of nonzero
%set smaller testnum to save time
testnum = 50; % number of tests
thresh = 1e-3;
lam = 1e-6; q = 0.5;
opts.s0 = round(m/2);
opts.gamma = 0.9;
maxit = 500; % max number of iterations
tol = 1e-4; % tolerance
len = length(k_seq);
succ = zeros(1,len);
for i = 1:len
    k = k_seq(i);
    for num = 1:testnum
        % generate a random test
        A = randn(m,n);
        x0 = zeros(n,1); 
        x0(randsample(n,k)) = randn(k,1);
        b = A*x0;        
        
        tic;
        [x,Out] = IRucLq_v(A,b,lam,q,opts); % call the vector-recovery solver
        time = toc; 
        
        relerr = norm(x-x0)/norm(x0);
        if relerr<=thresh 
            succ(i) = succ(i)+1;
        end
        fprintf('relerr = %4.3e, time = %4.3e\n', relerr, time);
    end
end
succ = succ/testnum; % get success rate
plot(k_seq,succ,'b-','linewidth',2);
xlabel('sparsity','fontsize',12);
ylabel('success rate', 'fontsize',12);
title('sparse vector recovery','fontsize',12)

%% low-rank matrix recovery
m = 100; n = 100; % matrix size
sr = 0.5; % sample ratio
p = round(sr*m*n); % number of samples
%set smaller testnum to save time
testnum = 20; % number of tests
lam = 1e-6; q = 0.5;
rank_set = 6:3:36;
len = length(rank_set);
succ_num = zeros(len,2);
thresh = 1e-3;
maxit = 1000; tol = 1e-5;

for rnum = 1:len
    r = rank_set(rnum);
    % over-estimate rank; of course, use the correct rank if you know it as it will improve the results
    esr = floor(1.5*r); % rank estimate, set to 150% of the true rank
    for num = 1:testnum
        M = randn(m,r)*randn(r,n);
        known = randsample(m*n,p);
        b = M(known);
        for sol_num = 1:2
            if (rnum>1 && succ_num(rnum-1,sol_num)>0) || rnum==1
                fprintf('r = %d, num = %d\n', r, num);
                opts = [];
                opts.maxit = maxit; opts.tol = tol;
                opts.rank = esr; opts.gamma = 0.9;
                opts.rank_adjust = 1; opts.min_rank = 5;
                
                t0 = tic;
                if sol_num==1
                    [X, Out] = IRucLq_m(m,n,known,b,lam,q,opts);
                else
                    [X, Out] = tIRucLq_m(m,n,known,b,lam,q,opts);
                end
                time = toc(t0);
                
                relres = norm(X-M,'fro')/norm(M,'fro');
                if relres<thresh
                    succ_num(rnum,sol_num) = succ_num(rnum,sol_num)+1;
                end
                fprintf('solver = %d: iter = %d, time = %4.2e, relres = %6.4e\n\n',...
                    sol_num, Out.iter, time, relres);
            end
        end
    end
    for sol_num = 1:2
        if succ_num(rnum,sol_num) == 0
            succ_num(rnum+1:end,sol_num) = 0;
        end
    end
end

succ_rate = succ_num'/testnum;
fig = figure('PaperSize',[5,4],'PaperPosition',[0,0,5,4]);
hold on;
plot(rank_set,succ_rate(1,:),'r--o');
plot(rank_set,succ_rate(2,:),'b-v');
axis([6,36,-0.1,1.1]);
legend('IRucLq-M','t-IRucLq-M','location','best');
xlabel('Rank','fontsize',12);
ylabel('Frequency of success','fontsize',12);
title('low-rank matrix recovery','fontsize',12);
