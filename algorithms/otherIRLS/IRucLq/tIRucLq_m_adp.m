function [Xr, outs] = tIRucLq_m(prob,lambda,opts)

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
% =========================================================================
% Adapted by Christian Kuemmerle, 
% Johns Hopkins University, kuemmerle@jhu.edu, 2020.
% Changes to original code by Yangyang Xu:
% - changed input/output interfaces.
% - added different update rules for smoothing parameter gam that uses the 
%   estimate of the model order (if opts.mode_eps == 'oracle_model_order')
%   or the relative change of two consecutive iterates (if opts.mode_eps 
%   == 'iter_diff').
% - parameter 'lambda' in this function corresponds to 'lam*q' in the
%   original code by Y. Xu. This adjustment was done to allow for a
%   generalization of the algorithm to q=0, and to allow for better
%   comparison with MatrixIRLS.
% =========================================================================
[r,rank_adjust,min_rank,Phi,X0_revealed,N0,saveiterates,verbose,...
    mode_eps,gamma0,gamma,q,tol]=get_options_problem_paras(prob,opts);
[d1,d2] = size(X0_revealed);

B = X0_revealed;
A = Phi;
I = speye(d1);
p = q/2-1; 
lamq = lambda; %lambda*q; (previous choice)
min_iter = 10; 
nstall = 0;
linopt.SYM = true; 
linopt.POSDEF = true;

E = ones(d1,d2); e = lamq*ones(d1,1);
loc = (1:(d1+1):d1^2)';

X_c = full(B./(lamq*E+A));

% X0 = X+ones(m,n);
eigopts.issym = 1;
eigopts.issreal = 1;
eigopts.disp = 0;

[V,D] = eigs(X_c*X_c',r+1,'lm',eigopts); 
d = diag(D);
eps_c = sqrt(d(r+1));
V = V(:,1:r); d = d(1:r); d = (d+eps_c^2).^p-eps_c^(q-2);

hist_d = zeros(1,N0);
dK = 1; dK0 = 1;

time = zeros(1,N0);
eps   = zeros(1,N0);
if saveiterates 
    X    = cell(1,N0);
    sings = cell(1,N0);
end
if verbose > 0
    fprintf(['Run IRucLq with q=%.2f, data fit parameter lambda=%.2e and smoothing parameter rule "%-11s"...\n'], ...
    q,lambda,mode_eps);
    dold = zeros(d1,1);
end
tic
%main iteration
for k = 1:N0
    eps(k) = eps_c;
    crit1 = dK<tol;
    crit2 = abs(dK-dK0)<tol*max(1,dK0);
    if crit2 
        nstall = nstall+1; 
    else
        nstall = 0;
    end
    if crit1 || nstall>=3
        switch crit1
            case 1
                 fprintf('\nThe %dth singular value is %f.\t Converged!\n',...
                     r+1, eps_c/gamma);
            case 0
                fprintf('\nNo more improvement\n')
        end
        if rank_adjust == 2 & r<K0
            [U,S,V] = svds(X_c,r);
            X_c = U*S*V';
        end
        break;
    end
    %updating X by solving equations
    D = spdiags(1./(lamq*d),0,r,r); 
    parfor ii=1:d2
        u = eps_c^(q-2)*e+A(:,ii); U = spdiags(1./u,0,d1,d1);         
        y = linsolve(-(D+V'*U*V),V'*(B(:,ii)./u),linopt);
        X_c(:,ii) = (B(:,ii)+V*y)./u;
    end    
    [V,D] = eigs(X_c*X_c',r+1,'lm',eigopts); 
    d = diag(D);
    if saveiterates
        X{k} = X_c;
        sings{k} = sort(sqrt(d(1:r)),'descend');
    end
    %adjust the rank estimate
    if rank_adjust == 1 && r>min_rank

        dR = d(min_rank:end-1);
        drops = dR(1:end-1)./dR(2:end);
        [dmx,imx] = max(drops);  
        rel_drp = (r-min_rank)*dmx/(sum(drops)-dmx);
        
        %If a large drop is found, we can adjust the rank estimate to imx 
        if rel_drp>10 
            fprintf('Rank adjust from %d to %d at iteration %d\n',...
                r, imx+min_rank-1, k);
            r = imx+min_rank-1; rank_adjust = 2;
        end
    end
    
    dK0 = dK;
    dK = sqrt(d(r+1));
    hist_d(k) = dK;
    if strcmp(mode_eps,'auto_decay')
        eps_c = max(gamma*eps_c,1e-12);
    elseif strcmp(mode_eps,'oracle_model_order')
        if k<=min_iter
            eps_c = min(eps_c,gamma0*sqrt(d(r+1)));
        else
            eps_c = min(eps_c,gamma*sqrt(d(r+1)));
        end
    elseif strcmp(mode_eps,'iter_diff')
        diff_c = norm(d(1:r)-dold(1:r));
        eps_c = min(eps_c,(diff_c/sqrt(d1))^2);
    end
    if verbose
        rel_chg = norm(d(1:r)-dold(1:r))./norm(d(1:r));
        fprintf(['k: %3d rankest: %3d eps: %.3e relchg: %.2e\n'], ...
        k,r,eps_c,rel_chg); 
        dold = d;
    end
    %use truncated eigen-decomposition
    V = V(:,1:r); 
    d = d(1:r);
    d = (d+eps_c^2).^p-eps_c^(q-2);
    time(k)=toc;
end
N = k-1;
fprintf('\n');
Xr = X_c;

outs = struct;
outs.N   = N;
outs.eps = eps(1:N);
outs.d   = hist_d(1:N);
outs.opts  = opts;
if saveiterates
    outs.X = X(1:N);
    outs.sings = sings(1:N);
end
outs.time=time(1:N);
end

function [r,rank_adjust,min_rank,Phi,X0_revealed,N0,saveiterates,verbose,...
    mode_eps,gamma0,gamma,q,tol]=get_options_problem_paras(prob,opts)
Phi = prob.Phi;
X0_revealed = prob.X0_revealed;

if isfield(opts,'rank_adjust')  
    rank_adjust = opts.rank_adjust;
else
    rank_adjust = 0;
end
if isfield(prob,'r') && not(isempty(prob.r)) && not(rank_adjust)
    r = prob.r;
    rank_adjust = 0;
    min_rank = r;
else
    r = [];
    %set rank_adjust 1 means to use a decreasing rank adjust strategy
    %min_rank is the minimum rank we can decrease the estimate rank to
    rank_adjust = 1; 
    if isfield(opts,'min_rank')     
        min_rank = opts.min_rank;
    else
        min_rank = 1;
    end
    r = round(min(d1,d2)/10);
end

% if isfield(opts,'W0')           
%     W = opts.W0;
% else
%     W = eye(d1);
% end
N0 = opts.N0;         
gamma0 = opts.gamma0;
gamma = opts.gamma;

if isfield(opts,'saveiterates') && opts.saveiterates == 1
    saveiterates = 1;
else
    saveiterates = 0;
end

mode_eps = opts.mode_eps;
q       = opts.p;   % non-convexity parameter
tol     = opts.tol;
verbose = opts.verbose;
end