function [x,Out] = IRucLq_v(A,b,lam,q,opts)

%Purpose: iterative reweighted method for unconstrained lq minimization for
%         sparse vector recovery
%Input:
%       A: sensing matrix
%       b: measurements
%       lam: balancing parameter between regularity and fidelity term
%       q: between (0,1]
%       opts.x0: initial point
%       opts.s0: estimated sparsity
%       opts.eps0: initial smoothing parameter
%       opts.tol: stopping tolerance
%       opts.maxit: maximum number of iterations
%       opts.gamma: control how fast eps is decreased
%Output:
%       x: recovered vector
%       Out.iter: number of iterations
%
%   copyright (c) by Yangyang Xu, Oct. 2012
%   yangyang.xu@rice.edu
%
% Paper:
% M.-J. Lai, Y. Xu, W. Yin. M.-J. Lai, Y. Xu, and W. Yin, Improved
% iteratively reweighted least squares for unconstrained smoothed lq
% minimization, SIAM Journal on Numerical Analysis, 5(2), 927-957, 2013.

[m,n] = size(A);
maxit = 500; tol = 1e-5;
eps = 1;  gamma = 0.9;
em = ones(1,m);  s = round(m/2);  x = zeros(n,1);

if isfield(opts,'x0')       x = opts.x0;                 end
if isfield(opts,'s0')       s = opts.s0;                 end
if isfield(opts,'eps0')     eps = opts.eps0;             end
if isfield(opts,'maxit')    maxit = opts.maxit;          end
if isfield(opts,'tol')      tol = opts.tol;              end
if isfield(opts,'gamma')    gamma = opts.gamma;          end

s0 = s; lamq = lam*q; p = 1-q/2;
At = A'; Atb = At*b; 

I = speye(m);
linopt.SYM = true; linopt.POSDEF = true;
nstall = 0; min_iter = 10; 

%denote the (s+1)th largest component
ds0 = 1; ds = 1;

for iter = 1:maxit
    
    %solve reweighted equation
    w = eps^2+x.^2; w = w.^p/lamq;
    wAtb = w.*Atb;
    wAt = (w*em).*At;
    w = linsolve(I+A*wAt, A*wAtb, linopt);
    x = wAtb-wAt*w;
    
    %get the (s+1)th largest component
    y = sort(abs(x),'descend');
    ds0 = ds;
    ds = y(s+1);
    
    %decrease smoothing parameter eps
    if iter<min_iter
        eps = min(eps,0.9*y(s+1));
    else eps = min(eps,gamma*y(s+1));
    end

    crit1 = ds<tol;
    crit2 = abs(ds-ds0)<tol*max(1,ds0);
    if crit2 nstall = nstall+1; else nstall = 0; end
    if crit1 || nstall>=3  break;        end
    
end
Out.iter = iter; 