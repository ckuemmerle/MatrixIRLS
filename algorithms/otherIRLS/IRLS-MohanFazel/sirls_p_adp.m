%% -------------- sIRLS-p (0 <= p <= 1) algorithm -------------------- %

%% ----- This is the code associated with the paper:
% ----- "Iterative Reweighted Algorithms for Matrix Rank Minimization"
% ----- Karthik Mohan (karna@uw.edu) and Maryam Fazel (mfazel@uw.edu).

% -------------- LAST UPDATE: 8/28/2012 ------------------------------ %
% =========================================================================
% Adapted by Christian Kuemmerle, 
% Johns Hopkins University, kuemmerle@jhu.edu, 2020.
% Changes to original code by Karthik Mohan and Maryam Fazel:
% - changed input/output interfaces.
% - added different update rules for smoothing parameter gam that uses the 
%   estimate of the model order (if opts.mode_eps == 'oracle_model_order')
%   or the relative change of two consecutive iterates (if opts.mode_eps 
%   == 'iter_diff')
% =========================================================================

function [Xr,outs] = sirls_p(prob,opts)
%% Get PARAMETERS
[r,rmax,rknown,Phi,X0_revealed,alpt,betat,N0,saveiterates,verbose,...
    mode_eps,incr,eta,gam0,gammin,p,tol]=get_options_problem_paras(prob,opts);
[d1,d2] = size(X0_revealed);
 
if(rknown == 1)
    countstart = r; %-- r if rank known
else
    countstart = rmax; %-- rmax if rank NOT known
end
normX0 = norm(full(X0_revealed));
if strcmp(mode_eps,'auto_decay') || strcmp(mode_eps,'iter_diff')
    gam = gam0*normX0;
elseif strcmp(mode_eps,'oracle_model_order')
    gam = normX0;
end

%% FIRST ITERATION OF sIRLS
k = 1;
Xnew = X0_revealed; 
Xold = X0_revealed;
svditer =  N0; %Parameter in rand_svd
count = countstart;
L = 2; %Initial Lipschitz constant
V = 0; D1 = 0;
extra_rank = 0;

time = zeros(1,N0);
eps   = zeros(1,N0);
gamma = zeros(1,N0);  
tic
if saveiterates 
    X    = cell(1,N0);%zeros(m,n,niter);
end
%% Run algorithm
if verbose > 0
    fprintf(['Run sIRLS-p with p=%.2f and smoothing parameter rule "%-11s"...\n'], ...
    p,mode_eps);
end
while(k<N0)
    [Xnew,~,terr,l] = grad_proj(X0_revealed,L,Xnew,V,D1,d1,d2,alpt,betat,2);
    if saveiterates
        X{k}=Xnew;
    end
    if strcmp(mode_eps,'auto_decay') || strcmp(mode_eps,'iter_diff')
        [U,S,V] = rand_svd(Xnew,count,k,svditer,incr);
        s = diag(S);
        V = V(:,1:size(s,1));
    elseif strcmp(mode_eps,'oracle_model_order')
        [U,S,V] = rand_svd(Xnew,count+1,k,svditer,incr);
        s = diag(S);
        s_long = s(r+1);
        s = s(1:count);
        V = V(:,1:count);
    end
    g = gam*ones(size(s,1),1); 
    s = s(1:size(s,1),1);
    D1 = diag( (g.^(1 - p/2))./((s.*s + g).^(1 - p/2))- ones(size(s,1),1));

    count = min(size(find(s > max(s)*1e-2),1)+extra_rank,rmax); 
    % Estimating rank to truncate SVD in next iteration
    if(rknown == 1)
        count = r;
    end

    L = 2; %Lipschitz constant
    diff_c = norm(Xnew - Xold,'fro');
    rel_chg = diff_c/norm(Xnew,'fro');
    Xold = Xnew;
    if strcmp(mode_eps,'auto_decay')
        gam = max(gam/eta,gammin);
    elseif strcmp(mode_eps,'oracle_model_order')
        gam = min(gam,(s_long/min(d1,d2)^(p/(2-p))).^2); %min(gam,(s_long/min(d1,d2)^q).^2);
    elseif strcmp(mode_eps,'iter_diff')
        gam = (diff_c/sqrt(min(d1,d2)))^2;
    end
    time(k) = toc;
    gamma(k) = gam;
    eps(k)   = sqrt(gam);
    if verbose
        eps_c = sqrt(gam);
        fprintf(['k: %3d gamma: %.3e relchg: %.2e\n'], ...
        k,gam,rel_chg);  
    else
        if(mod(k,20) == 0)
            fprintf('.');
        end
    end
    k = k + 1;
    if(mod(k,20) == 0)
        if(rel_chg < tol)
            break;
        end
    end
end
%% Prepare output
N = k-1;
Xr = Xnew;

outs = struct;
outs.N = N;
outs.eps = eps(1:N);
outs.gamma = gamma(1:N);
outs.opts  = opts;
if saveiterates
    outs.X = X(1:N);
end
outs.time=time(1:N);
end

function [r,rmax,rknown,Phi,X0_revealed,alpt,betat,N0,saveiterates,verbose,...
    mode_eps,incr,eta,gam0,gammin,p,tol]=get_options_problem_paras(prob,opts)
if isfield(prob,'r') && not(isempty(prob.r))
    rknown = 1;
    r = prob.r;
else
    rknown = 0;
    r = [];
end
if isfield(opts,'R')
    rmax = opts.R; % Upper bound for number of singular values to be calculated
else
    rmax = Inf;
end
Phi = prob.Phi;
X0_revealed = prob.X0_revealed;
[alpt,betat]= find(Phi);

N0 = opts.N0;


if isfield(opts,'saveiterates') && opts.saveiterates == 1
    saveiterates = 1;
else
    saveiterates = 0;
end

mode_eps = opts.mode_eps;
incr    = opts.incr;
eta     = opts.eta;
gam0    = opts.gam0;
gammin  = opts.gammin;
p       = opts.p;   % non-convexity parameter
tol     = opts.tol;
verbose = opts.verbose;
end
   