function [X, Out] = IRucLq_m(m,n,known,b,lam,q,opts)

%Purpose: iterative reweighted method for unconstrained lq minimization for
%         low-rank matrix recovery
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

%   copyright (c) by Yangyang Xu, Oct. 2012
%   yangyang.xu@rice.edu
%
% Paper:
% M.-J. Lai, Y. Xu, W. Yin.
% Improved iteratively reweighted least squares for unconstrained smoothed lq
% minimization, SIAM Journal on Numerical Analysis, 5(2), 927-957, 2013.

K = round(min(m,n)/10); W = eye(m); eps = 1;

gamma0 = 0.9;  gamma = 0.9;

maxit = 500; tol = 1e-5;

%set rank_adjust 1 means to use a decreasing rank adjust strategy
%min_rank is the minimum rank we can decrease the estimate rank to
rank_adjust = 0; min_rank = 1;

if isfield(opts,'rank')         K = opts.rank;                  end
if isfield(opts,'W0')           W = opts.W0;                    end
if isfield(opts,'eps0')         eps = opts.eps0;                end
if isfield(opts,'maxit')        maxit = opts.maxit;             end
if isfield(opts,'tol')          tol = opts.tol;                 end
if isfield(opts,'rank_adjust')  rank_adjust = opts.rank_adjust; end
if isfield(opts,'min_rank')     min_rank = opts.min_rank;       end
if isfield(opts,'gamma0')       gamma0 = opts.gamma0;           end
if isfield(opts,'gamma')        gamma = opts.gamma;             end



B = zeros(m,n); B(known) = b;
A = zeros(m,n); A(known) = 1; 
I = speye(m);

p = q/2-1; lamq = lam*q; min_iter = 10; 
nstall = 0;
linopt.SYM = true; linopt.POSDEF = true;

E = ones(m,n); 
Eps = zeros(1,maxit);
loc = (1:(m+1):m^2)';

X = B./(lamq*E+A);

[V,D] = eig(X*X'); d = diag(D);
eps0 = eps; eps = sqrt(d(m-K));
d = lamq*(d+eps^2).^(q/2-1);
W = V*diag(d)*V';

hist_d = zeros(1,maxit);
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
                     K+1, 2*eps);
            case 0
                fprintf('\nNo more improvement\n')
        end
        if rank_adjust == 2
            [U,S,V] = svds(X,K);
            X = U*S*V';
        end
        break;
    end
    
    %updating X by solving equations
    parfor j=1:n
        Cw = W; Cw(loc) = Cw(loc)+ A(:,j);
        R = chol(Cw);
        X(:,j) = R\(R'\B(:,j));
    end

    
    [V,D] = eig(X*X'); d = diag(D); 
    d = max(0,d); % for stability
    
    %adjust the rank estimate
    if rank_adjust == 1 & K>min_rank

        dR = d(end-min_rank+1:-1:end-K+1);
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
    dK = sqrt(d(m-K));
    hist_d(iter) = sqrt(d(m-K));
    eps0 = eps;
    if iter<=min_iter
        eps = min(eps,gamma0*sqrt(d(m-K)));
    else
        eps = min(eps,gamma*sqrt(d(m-K)));
    end
    
    d = lamq*(d+eps^2).^p;
    W = V*diag(d)*V';
end
fprintf('\n');
Eps = Eps(1:iter);
Out.iter = iter; Out.eps = Eps; Out.d = hist_d(1:iter-1);
end
