function [Xr,outs] = ScaledGD(prob,opts,r)
%
% Implementation of Scaled Gradient Descent ('ScaledGD') for
% matrix completion [1].
% Code by Tian Tong, cf. https://github.com/Titan-Tong/ScaledGD.
% =========================================================================
% [1] Tian Tong, Cong Ma, Yuejie Chi. Accelerating Ill-Conditioned 
% Low-Rank Matrix Estimation via Scaled Gradient Descent 
% overrelaxation algorithm, preprint, arXiv:2005.08898.
% =========================================================================
% This basically corresponds to the implementation of the algorithm
% provided in "MC.m" by Tong et al.
% Modifications by Christian Kuemmerle:
% - save intermediate iterates if opts.saveiterates == 1.
% - save timing information.
% - increase normalization factor 'p' by a factor of 1.5 each time the
% algorithmic iterates diverge for a given p.

d1=prob.d1;
d2=prob.d2;
% if not(d1 == d2)
%    error('Not implemented yet.') 
% end
m = length(prob.Omega);
T = opts.T; % number of gradient steps
 
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    saveiterates = 1;
    Xouts = cell(1,T);
else
    saveiterates = 0;
end
if isfield(opts,'scaledgd_normalization')
    p = opts.scaledgd_normalization;
else
    p = m/(d1*d2); % normalization parameter
end
T = opts.T; % number of gradient steps
eta = opts.eta; % step size parameter

%% Initialization
% U_seed = sign(rand(d1,r)-0.5);
% [U_star,~,~] = svds(U_seed, r);
% V_seed = sign(rand(d2,r)-0.5);  
% [V_star,~,~] = svds(V_seed, r);
% d = min(d1,d2);

%mu = max(vecnorm([U_star;V_star], 2, 2))^2*d/r; 


p_good = 0;
p=p/1.5;
while p_good == 0
    p_good = 1;
    p=p*1.5;
    X_sps=sparse(prob.rowind,prob.colind,ones(1,m),d1,d2);
    Y_mat=sparse(prob.rowind,prob.colind,prob.y,d1,d2);
    [U0, Sigma0, V0] = svds(Y_mat/p, r);
    L = U0*sqrt(Sigma0);
    R = V0*sqrt(Sigma0);
    X_vec = partXY(L',R',prob.rowind,prob.colind,m);
    setSval(X_sps,X_vec,m);
    norm_X0Omega=norm(prob.y,2);
    outs = struct;
    time = zeros(1,T);
    if saveiterates
        Xouts{1} = {L,R};
    end
    tic

    %% Run Scaled Gradient Descent
    for t = 1:T
        %X = L*R';
    %     Xouts{t} = {L,R};
    % %     dist_X = norm(X - X_star, 'fro')/norm(X_star, 'fro');
    % %     dists_X_scaledGD(i_kappa, t, trial) = dist_X;
    %     if dist_X < 1e-14
    %         break;
    %     end
        Z = (X_sps - Y_mat)/p; %(Omega.*X - Y_mat)/p;
        L_plus = L - eta*Z*R/(R'*R + eps('double')*eye(r));
        R_plus = R - eta*Z'*L/(L'*L + eps('double')*eye(r));
        
        if norm(L_plus)>2*norm(L) || norm(R_plus)>2*norm(R)
            %norm_fac_L=norm(L_plus)/norm(L);
            %norm_fac_R=norm(R_plus)/norm(R);
            p_good = 0;
            if opts.verbose > 0
                disp('Update p factor...')
            end
            break;
        end
        L = L_plus;
        R = R_plus;
        X_vec = partXY(L',R',prob.rowind,prob.colind,m);
        setSval(X_sps,X_vec,m);
        if saveiterates
            Xouts{t} = {L,R};
        end
        time(t) = toc;
        if t > 1 
            rel_dist = norm(X_vec-prob.y.',2)/norm_X0Omega;
            if rel_dist < opts.tol
                break;
            end
        end
    end
end
outs.time=time(1:t);
outs.p = p;
if saveiterates
    Xr = Xouts(1:t);
else
    Xr = {{L,R}};
end
end

