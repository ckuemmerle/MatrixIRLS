clear, close all, clc;
n = 1000; 
r = 10;
kappa_list = [1,5,10,20];
alpha = 0.1;

T = 1000;
eta = 0.5;
len_trial = 1;
dists_X_scaledGD = zeros(length(kappa_list), T, len_trial);
dists_X_GD = zeros(length(kappa_list), T, len_trial);

for trial = 1:len_trial
    U_seed = sign(rand(n,r)-0.5);
    [U_star,~,~] = svds(U_seed, r);
    V_seed = sign(rand(n,r)-0.5);  
    [V_star,~,~] = svds(V_seed, r);
    mu = max(vecnorm([U_star;V_star], 2, 2))^2*n/r; 
    S_seed = randn(n, n);
    S_star = Tproj(S_seed, alpha);    
    for i_kappa = 1:length(kappa_list)    
        kappa = kappa_list(i_kappa);
        sigma_star = linspace(kappa, 1, r);
        L_star = U_star*diag(sqrt(sigma_star));
        R_star = V_star*diag(sqrt(sigma_star));
        X_star = L_star*R_star';
        Y = X_star + S_star;
        [U0, Sigma0, V0] = svds(Y - Tproj(Y, alpha), r);
        %% Scaled GD
        L = U0*sqrt(Sigma0);
        R = V0*sqrt(Sigma0);
        for t = 1:T
            X = L*R';
            dist_X = norm(X - X_star, 'fro')/norm(X_star, 'fro');
            dists_X_scaledGD(i_kappa, t, trial) = dist_X;
            if dist_X < 1e-14
                break;
            end
            S = Tproj(Y - X, 2*alpha);
            L_plus = L - eta*(X+S-Y)*R/(R'*R + eps('double')*eye(r));
            R_plus = R - eta*(X+S-Y)'*L/(L'*L + eps('double')*eye(r));
            L = L_plus;
            R = R_plus;
        end
        %% GD
        L = U0*sqrt(Sigma0);
        R = V0*sqrt(Sigma0);
        for t = 1:T
            X = L*R';
            dist_X = norm(X - X_star, 'fro')/norm(X_star, 'fro');
            dists_X_GD(i_kappa, t, trial) = dist_X;
            if dist_X < 1e-14
                break;
            end
            S = Tproj(Y - X, 2*alpha);
            L_plus = L - eta/sigma_star(1)*(X+S-Y)*R;
            R_plus = R - eta/sigma_star(1)*(X+S-Y)'*L;
            L = L_plus;
            R = R_plus;
        end
    end
end

figure('position', [200,200,800,800]);
clrs = {[.5,0,.5], [1,.5,0], [1,0,0], [0,1,0], [0,0,1]};
mks = {'o', 'x', 'p', 's', 'd'};
lgds = {};
for i_kappa = 1:length(kappa_list)
    kappa = kappa_list(i_kappa);
    dists = mean(dists_X_scaledGD(i_kappa, :, :), 3);
    dists = dists(dists > 1e-14);
    T_subs = (2*i_kappa):10:length(dists);
    semilogy(T_subs, dists(T_subs), 'Color', clrs{1}, 'Marker', mks{i_kappa}, 'MarkerSize', 9);
    hold on; grid on;
    lgds{end+1} = sprintf('$\\mathrm{ScaledGD}~\\kappa=%d$', kappa);
end
for i_kappa = 1:length(kappa_list)
    kappa = kappa_list(i_kappa);
    dists = mean(dists_X_GD(i_kappa, :, :), 3);
    dists = dists(dists > 1e-14);
    T_subs = 10:10:length(dists);
    semilogy(T_subs, dists(T_subs), 'Color', clrs{2}, 'Marker', mks{i_kappa}, 'MarkerSize', 9);
    hold on; grid on;
    lgds{end+1} = sprintf('$\\mathrm{VanillaGD}~\\kappa=%d$', kappa);
end
xlabel('Iteration count');
ylabel('Relative error');
legend(lgds, 'Location', 'Best', 'Interpreter', 'latex');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
fig_name = sprintf('RPCA_n=%d_r=%d_alpha=%g',n,r,alpha);

function S = Tproj(S, alpha)
% S=T_alpha(S) operator for sparse projection
    [n1, n2] = size(S);
    krow = floor(alpha*n1);
    kcol = floor(alpha*n2);
    [~, loc_row] = maxk(abs(S'), krow, 1);
    [~, loc_col] = maxk(abs(S), kcol, 1);
    id_row = repmat((1:n1), krow, 1);
    id_col = repmat((1:n2), kcol, 1);
    mask_row = sparse(id_row(:), loc_row(:), 1, n1, n2, nnz(S));
    mask_col = sparse(loc_col(:), id_col(:), 1, n1, n2, nnz(S));    
    S = S.*mask_row.*mask_col;
end