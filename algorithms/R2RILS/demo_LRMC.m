% simple demo for low rank matrix completion

iter = 1; 

m = 1000; n = m;%m+100; 

rho = 1.5; %oversampling ratio

r = 10; sigma_list = [1 1 1 1 1 1 1 1 1 1]; 
modeX0 = 2;
complex = 0;
nv = floor(r*(n+m-r) * rho);  %number of observed entries
t_max = 40;
opts.N0=t_max;
opts.assemblesparseA = 0;

for counter =1:iter
    % generate low rank matrix X with singular values in sigma_list 
    %[X0 Utrue Vtrue] = generate_low_rank_matrix(m,n,sigma_list);  %% Missing in Bauch's and Nadler's code
    [Utrue,Vtrue] = sample_X0_lowrank(m,n,r,modeX0,complex);
    X0 = Utrue*Vtrue';
    % generate random set of observed entries
    t = randperm(n*m); 
    omega = zeros(nv,2); 

    size_matrix = [m n]; 
    [omega(:,1) omega(:,2)] = ind2sub(size_matrix,t(1:nv)); 
    
    % generate matrix X with entries equal to X0(i,j) for (i,j) in omega
    X = zeros(m,n); 
    for i=1:nv
        X(omega(i,1),omega(i,2)) = X0(omega(i,1),omega(i,2));
    end
    tic; 
    [X_hat outs, observed_RMSE] = R2RILS(X,omega,r,opts); 
    elapsed_time = toc; 
    
    fprintf('iter %4d RMSE %8d\n',counter,sqrt( sum(sum((X_hat - X0).^2)) ) / sqrt(n*m) ); 
    fprintf('TIME %5.1f\n',elapsed_time); 
end
    
