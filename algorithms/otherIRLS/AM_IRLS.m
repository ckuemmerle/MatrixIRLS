function [X,N,eps,singval,time] = ...
    AM_IRLS(prob,opts)
%AM-IRLS Arithmetic mean iteratively least squares algorithm, implemented for
%matrix completion measurements.
%
%   Inputs: prob = Struct with fields "d1", "d2", "r", "Omega", "Phi" for
%                       problem speficiation
%           opts = Struct with fields "N0", "p", "epsmin", "tol" and
%                       "linsolve"
%  Outputs:    X = (1 x N) cell with (d1 x d2) matrices: 
%                       iterates of HM-IRLS algorithm
%              N = number of performed iterations
%            eps = (1 x N) vector with values of epsilon-smoothing at every
%                   iteration
%        singval = (1 x N) cell with vectors of singular values of X{l}
%           time = elapsed time
%
% For more information, see the paper:
% [1] C. Kuemmerle, J. Sigl, "Harmonic Mean Iteratively Reweighted Least
% Squares for Low Rank Matrix Recovery", Journal of Machine Learning 
% Research, 19(47):1-49, 2018.

d1          = prob.d1;       % first dimension of matrix to be recovered
d2          = prob.d2;       % second dimension of matrix to be recovered
r           = prob.r;        % target rank
Omega       = prob.Omega;    % indicates "seen" linear indices of matrix
Phi         = prob.Phi;      % measurement mask
y           = prob.y;        % measured entries

N0          = opts.N0;      % max. number of iterations
p           = opts.p;       % Schatten-p parameter
epsmin      = opts.epsmin;  % minimal eps
tol         = opts.tol;     % stopping criterion (minimal relative change of Frobenius norm)

dd=min(d1,d2);

%%% get row and column indices of entry mask
[rowind_mask,colind_mask]   = find(Phi);

%% Iterations
X                = cell(1,N0);
eps              = zeros(1,N0);
singval          = cell(1,N0);
time             = zeros(1,N0);

l=1;
eps_c = Inf;
U_X_c = eye(d1);
V_X_c = eye(d2);
H_c   = ones(d1,d2);
X_c   = zeros(d1,d2);
tic
while l <= N0
    X_cm1=X_c;
    if mod(l,50) == 0 
        disp(['Begin Iteration ',num2str(l),'...']);
    end
    
    if strcmp(opts.linsolve,'direct_solve')
        %% Update matrix of linear system to be solved
        PhWiPh_star = update_PhWiPh(H_c,U_X_c,V_X_c,Omega,rowind_mask,colind_mask);
        %% Calculate new iterate
            X_c = directsolve_weightedLS(PhWiPh_star,y,...
               U_X_c,V_X_c,Omega,'arithmetic',H_c);
    end
    %% Calculate new weights
    [U_X_c,singval_mat_c,V_X_c] = svd(X_c,'econ');
    singval_c = diag(singval_mat_c);
    eps_c = min(eps_c,max(singval_c(r+1),epsmin));
    Sigma_c_diag= (singval_c.^2+ eps_c.^2.*ones(dd,1)).^((2-p)/2);
    
    H_c=(kron(Sigma_c_diag(1:d1).^(-1),0.5.*ones(1,d2))+kron(0.5.*ones(d1,1),Sigma_c_diag(1:d2).^(-1)')).^(-1);
    eps(l)      = eps_c;
    singval{l}  = singval_c;
    X{l}        = X_c;
    %% Check stopping criterion
   
    if l >= 2 && l <= N0 && norm(X_c-X_cm1,'fro')/norm(X_c,'fro') < tol
        N0 = l;
    end
    time(l)=toc;
    l=l+1;
end
%% Tidy up the size of the relevant output arrays and cells
X             = X(1:N0);
eps           = eps(1,1:N0);
singval       = singval(1,1:N0);
N=N0;
time          = time(1:N0);
end



function PhWiPh_star = update_PhWiPh(H_c,U_c,V_c,Omega,rowind_mask,colind_mask)

m =length(rowind_mask);
PhWiPh_star = zeros(m,m);
for ll = 1:m
    M_c=H_c.*(U_c(rowind_mask(ll),:)'*V_c(colind_mask(ll),:));
    M_c=(U_c*M_c)*V_c';
    PhWiPh_star(:,ll)=M_c(Omega);
end

end