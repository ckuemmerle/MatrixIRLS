function [X, Out] = tIRucLq_m(m,n,known,b,lam,q,opts)

%Purpose: iterative reweighted method for unconstrained lq minimization for
%         low-rank matrix recovery
%
%This code differs from IRucLq_m.m by applying a truncation acceleration;
%see the paragraph around (4.2) in our paper
%
%Input:
%       [m,n]: size of matrix
%       b: measurements
%       lam: balancing parameter between regularity and fidelity term
%       q: between (0,1]
%       opts.rank: estimated rank
%       opts.eps0: initial smoothing parameter
%       opts.tol: stopping tolerance
%       opts.maxit: maximum number of iterations
%       opts.gamma0: control how fast eps is decreased before min_iter
%       opts.gamma: control how fast eps is decreased after min_iter
%                   set a smaller gamma can accelerate convergence
%                   but be careful to choose gamma since the problem is
%                   nonconvex as q<1
%Output:
%       X: recovered matrix
%       Out.iter: number of iterations
%       Out.eps: history of the changes of smoothing parameter eps
%       Out.d: history of the changes of (K+1)th singular values
%
%   copyright (c) by Yangyang Xu, Oct. 2012
%   yangyang.xu@rice.edu
%
% Paper:
% M.-J. Lai, Y. Xu, W. Yin.
% Improved iteratively reweighted least squares for unconstrained smoothed lq
% minimization, SIAM Journal on Numerical Analysis, 5(2), 927-957, 2013. 

K = round(min(m,n)/10); W = eye(m); eps = 1;  

gamma0 = 0.9; gamma = 0.9;

maxit = 500; tol = 1e-5;

%set rank_adjust 1 means to use a decreasing rank adjust strategy
%min_rank is the minimum rank we can decrease the estimate rank to
rank_adjust = 1; min_rank = 1;

if isfield(opts,'rank')        K = opts.rank;                   end
if isfield(opts,'W0')          W = opts.W0;                     end
if isfield(opts,'eps0')        eps = opts.eps0;                 end
if isfield(opts,'maxit')       maxit = opts.maxit;              end
if isfield(opts,'tol')         tol = opts.tol;                  end
if isfield(opts,'rank_adjust') rank_adjust = opts.rank_adjust;  end
if isfield(opts,'min_rank')    min_rank = opts.min_rank;        end
if isfield(opts,'gamma0')      gamma0 = opts.gamma0;            end
if isfield(opts,'gamma')       gamma = opts.gamma;              end

K0 = K;

B = zeros(m,n); B(known) = b;
A = zeros(m,n); A(known) = 1; 

p = q/2-1; lamq = lam*q; min_iter = 10; 
nstall = 0;
linopt.SYM = true; linopt.POSDEF = true;

E = ones(m,n); e = lamq*ones(m,1);

X = B./(lamq*E+A); X0 = X+ones(m,n);
eigopts.issym = 1;
eigopts.issreal = 1;
eigopts.disp = 0;

[V,D] = eigs(X*X',K+1,'lm',eigopts); d = diag(D);
eps0 = eps; eps = sqrt(d(K+1));
V = V(:,1:K); d = d(1:K); d = (d+eps^2).^p-eps^(q-2);

Eps = zeros(1,maxit);

dK = 1; dK0 = 1;

%main iteration
for iter = 1:maxit
    
    Eps(iter) = eps;
    
    crit1 = dK<tol;
    crit2 = abs(dK-dK0)<tol*max(1,dK0);
    if crit2 nstall = nstall+1; else nstall = 0; end
    if crit1 || nstall>=3
        switch crit1
            case 1
                 fprintf('\nThe %dth singular value is %f.\t Converged!\n',...
                     K+1, eps/gamma);
            case 0
                fprintf('\nNo more improvement\n')
        end
        if rank_adjust == 2 & K<K0
            [U,S,V] = svds(X,K);
            X = U*S*V';
        end
        break;
    end
    %updating X by solving equations
    D = spdiags(1./(lamq*d),0,K,K); 
    parfor ii=1:n
        u = eps^(q-2)*e+A(:,ii); U = spdiags(1./u,0,m,m);         
        y = linsolve(-(D+V'*U*V),V'*(B(:,ii)./u),linopt);
        X(:,ii) = (B(:,ii)+V*y)./u;
    end    
    
    [V,D] = eigs(X*X',K+1,'lm',eigopts); d = diag(D);
    
    %adjust the rank estimate
    if rank_adjust == 1 & K>min_rank

        dR = d(min_rank:end-1);
        drops = dR(1:end-1)./dR(2:end);
        [dmx,imx] = max(drops);  
        rel_drp = (K-min_rank)*dmx/(sum(drops)-dmx);
        
        %If a large drop is found, we can adjust the rank estimate to imx 
        if rel_drp>10 
            fprintf('Rank adjust from %d to %d at iteration %d\n',...
                K, imx+min_rank-1, iter);
            K = imx+min_rank-1; rank_adjust = 2;
        end
    end
    
    dK0 = dK;
    dK = sqrt(d(K+1));
    
    if iter<=min_iter
        eps = min(eps,gamma0*sqrt(d(K+1)));
    else
        eps = min(eps,gamma*sqrt(d(K+1)));
    end
    
    %use truncated eigen-decomposition
    V = V(:,1:K); d = d(1:K); d = (d+eps^2).^p-eps^(q-2);
end
fprintf('\n');
Eps = Eps(1:iter);
Out.iter = iter; Out.eps = Eps;
end
