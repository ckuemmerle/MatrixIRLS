function [X_hat U_hat lambda_hat V_hat, observed_RMSE] = R2RILS(X,omega,r,t_max)
%
% Implementation of Rank 2r iterative least squares ('R2RILS') for
% matrix completion [1].
% Code by Jonathan Bauch and Boaz Nadler, 2020, cf.
% - https://github.com/Jonathan-WIS/R2RILS
% - http://www.wisdom.weizmann.ac.il/~nadler/Projects/R2RILS/R2RILS.html
% =========================================================================
% [1] Jonathan Bauch and Boaz Nadler. Rank 2r iterative least squares: 
% efficient recovery of ill-conditioned low rank matrices from few entries, 
% preprint, arXiv:2002.01849.
% =========================================================================
% WRITTEN BY BAUCH & NADLER / 2020
%
% INPUT: 
% X = matrix with observed entries in the set omega
% omega = list of size nv * 2 of pairs (i,j) ; nv=total number of visible entries
% r = target rank of reconstructed matrix
% t_max = maximal number of iterations [optional] 
%
%


eps_early_stop = 1e-15; 

if nargin < 4
    t_max = 50; % default total number of iterations
end

[nr nc] = size(X);   %nr,nc = number of rows / colums

nv = size(omega,1); 

Xmax = max(max(abs(X)));   %max absolute value of all observed entries

rhs = zeros(nv,1);   %vector of visible entries in matrix X

if 0 fprintf('Inside R2RILS nr= %d nc = %d nv= %d\n',nr,nc,nv);  end

for counter=1:nv
    rhs(counter) = X(omega(counter,1),omega(counter,2)); 
end

% Initialization by rank-r SVD of observed matrix
[U S V] = svds(X,r); % U is of size nr x rank; V is of size nc x rank (both column vectors)

U0 = U; V0 = V; 

m = (nc + nr) * r;  % total number of variables in single update iteration

% Z^T = [a(coordinate 1)  a(coordinate 2) ... a(coordinate nc) | b(coordinate 1) ... b(coordinate nr) ]
% Z^T = [ Vnew and then Unew ] 

observed_RMSE=zeros(t_max,1); 

X_hat_previous = zeros(nr,nc); 

for loop_idx = 1:t_max
    
    if 0 fprintf('loop_idx  %d/%d\n',loop_idx,t_max); end

    Z = zeros(m,1);    % Z is a long vector with a vectorization of A and of B (notation in paper)
    A = zeros(nv,m);   % matrix of least squares problem 

    %CONSTRUCTION OF MATRIX A, 
    for counter=1:nv
        j = omega(counter,1); k = omega(counter,2); 
        index=r*(k-1)+1; 
        A(counter,index:index+r-1) = U(j,:); 
        index = r*nc + r*(j-1)+1; 
        A(counter,index:index+r-1) = V(k,:); 
    end

    A = sparse(A); 
    
    Z = lsqr(A,rhs,1e-15,150);   % LSQR finds the minimum norm solution and is much faster than lsqminnorm
    
    % construct U and V from the entries of the long vector Z 
    Utilde = zeros(size(U)); Vtilde = zeros(size(V)); 
    nc_list = r* [0:1:(nc-1)]; 
    for i=1:r
        Vtilde(:,i) = Z(i+nc_list); 
    end

    nr_list = r*[0:1:(nr-1)]; 
    start_idx = r*nc; 
    for i=1:r
        Utilde(:,i) = Z(start_idx + i + nr_list);
    end
   
    X_hat = U *Vtilde' + Utilde * V';   % rank 2*r intermediate result

    observed_RMSE(loop_idx) = sqrt(sum(sum( (abs(X)>0).*(X_hat-X).^2) )  / nv); 
    
    normU = sqrt(sum(Utilde.^2)); 
    Utilde = Utilde * diag(1./normU);  

    normV = sqrt(sum(Vtilde.^2)); 
    Vtilde = Vtilde * diag(1./normV);  

    % AVERAGING STEP FOLLOWED BY COLUMN NORMALIZATION
    U = (U + Utilde);
    V = (V + Vtilde); 
    
    normU = sqrt(sum(U.^2)); 
    U = U * diag(1./normU); 

    normV = sqrt(sum(V.^2)); 
    V = V * diag(1./normV); 

    fprintf('INSIDE R2RILS %3d DIFF X_hat %8d\n',loop_idx,sqrt(sum(sum((X_hat-X_hat_previous).^2)))/sqrt(nr*nc)); 
    %EARLY STOPPING
    if sqrt(sum(sum((X_hat-X_hat_previous).^2)))/sqrt(nr*nc) <eps_early_stop
        break; 
    end
    
    X_hat_previous = X_hat ;
end

% return rank r SVD of X_hat which originally is of rank 2r
[U_hat lambda_hat V_hat] = svds(X_hat,r); 

X_hat = U_hat * lambda_hat * V_hat'; 
