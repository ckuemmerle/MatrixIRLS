function [X_hat, outs, observed_RMSE] = R2RILS(X,omega,r,opts)
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
% % Modifications by Christian Kuemmerle:
% - save intermediate iterates if opts.saveiterates == 1.
% - save timing information.
% - improved assembling of linear system matrix A using sparse matrix
% construction.
% - more flexible algoritmic options (see opts).
%
% INPUT: 
% X = matrix with observed entries in the set omega
% omega = list of size nv * 2 of pairs (i,j) ; nv=total number of visible entries
% r = target rank of reconstructed matrix
% opts = struct with algorithmic parameters. Field may include:
%       - N0_inner: Max. nr. of iterations
%       - tol_CG_fac: Tolerance for stopping condition of iterative linear
%       system solver (LSQR)
%       - tol: Tolerance for stopping condition of outer iterations
%
verbose = opts.verbose;

%eps_early_stop = 1e-15;  choice by authors
eps_early_stop = opts.tol;
tol_inner      = 1e-1*opts.tol_CG_fac;
if isfield(opts,'N0_inner')
    N0_inner = opts.N0_inner;
else
    N0_inner = 150;
end
t_max = opts.N0;

% if nargin < 4
%     t_max = 50; % default total number of iterations
% end
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    Xout=cell(1,t_max);% Xout=zeros(prob.n1,prob.n2,opts.maxit);
end
assemblesparseA = opts.assemblesparseA;

time=zeros(1,t_max);

[nr nc] = size(X);   %nr,nc = number of rows / colums

nv = size(omega,1); 

Xmax = max(max(abs(X)));   %max absolute value of all observed entries

rhs = zeros(nv,1);   %vector of visible entries in matrix X

if verbose > 0
    if 0 fprintf('Inside R2RILS nr= %d nc = %d nv= %d\n',nr,nc,nv);  end
end

for counter=1:nv
    rhs(counter) = X(omega(counter,1),omega(counter,2)); 
end
if assemblesparseA
   colind_A = generate_sparse_matrix_indices(omega,r,nr,nc); 
end
% Initialization by rank-r SVD of observed matrix
[U S V] = svds(X,r); % U is of size nr x rank; V is of size nc x rank (both column vectors)

U0 = U; V0 = V; 

m = (nc + nr) * r;  % total number of variables in single update iteration

% Z^T = [a(coordinate 1)  a(coordinate 2) ... a(coordinate nc) | b(coordinate 1) ... b(coordinate nr) ]
% Z^T = [ Vnew and then Unew ] 

observed_RMSE=zeros(t_max,1); 

X_hat_previous = zeros(nr,nc); 

tic
for loop_idx = 1:t_max
    if verbose > 0
        if 0 fprintf('loop_idx  %d/%d\n',loop_idx,t_max); end
    end

    Z = zeros(m,1);    % Z is a long vector with a vectorization of A and of B (notation in paper)
    if assemblesparseA
        A = generate_sparse_A(U,V,omega,colind_A);
    else
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
    end 
    if verbose > 0
        Z = lsqr(A,rhs,tol_inner,N0_inner);   % LSQR finds the minimum norm solution and is much faster than lsqminnorm
    else
        [Z,flag] = lsqr(A,rhs,tol_inner,N0_inner); 
    end
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
    rel_chg = sqrt(sum(sum((X_hat-X_hat_previous).^2))/sum(sum(X_hat.^2)));
    if verbose
        fprintf('INSIDE R2RILS %3d DIFF X_hat %8d\n',loop_idx,rel_chg);
    end
    if isfield(opts,'saveiterates') && opts.saveiterates == 1
        Xout{loop_idx}= X_hat;
    end
    time(loop_idx)=toc;
    %EARLY STOPPING
    if rel_chg < eps_early_stop
        break; 
    end
    
    X_hat_previous = X_hat ;
end
outs = struct;
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    outs.Xout = Xout(1:loop_idx);
end
outs.N=loop_idx;
outs.time = time(1:loop_idx);
% return rank r SVD of X_hat which originally is of rank 2r

[U_hat lambda_hat V_hat] = svds(X_hat,r); 
X_hat = U_hat * lambda_hat * V_hat'; 

end


function A = generate_sparse_A(U,V,omega,colind_A)
%%%
% Replicates: CONSTRUCTION OF MATRIX A
%     for counter=1:nv
%         j = omega(counter,1); k = omega(counter,2); 
%         index=r*(k-1)+1; 
%         A(counter,index:index+r-1) = U(j,:); 
%         index = r*nc + r*(j-1)+1; 
%         A(counter,index:index+r-1) = V(k,:); 
%     end
%%%
nv = size(omega,1); 
r = size(U,2);
nr = size(U,1);
nc = size(V,1);
%A = zeros(nv,m);   % matrix of least squares problem 
val_A=zeros(1,2*r*nv);
rowind_A=kron(1:nv,ones(1,2*r));
for counter=1:nv
    j = omega(counter,1); k = omega(counter,2);
    %colind_A = [colind_A,(r*(k-1)+1):(r*(k-1)+r),(r*nc + r*(j-1)+1):(r*nc + r*(j-1)+r)];
    val_A((counter-1)*(2*r)+1:counter*2*r) = [U(j,:),V(k,:)];
end
A = sparse(rowind_A,colind_A,val_A,nv,r*(nr+nc));
end

function colind_A = generate_sparse_matrix_indices(omega,r,nr,nc)

nv = length(omega);
colind_A=zeros(1,2*r*nv);
for counter=1:nv
    j = omega(counter,1); k = omega(counter,2);
    colind_A((counter-1)*(2*r)+1:counter*2*r) = [(r*(k-1)+1):(r*(k-1)+r),...
        (r*nc + r*(j-1)+1):(r*nc + r*(j-1)+r)];
end
end