clear, close all, clc;
n1 = 1000; n2=1000; nd=n1+n2-1; 
r = 10;
kappa_list = [1,5,10,20];
p = 0.2;

T = 1000;
eta = 0.5;
len_trial = 1;
dists_X_scaledGD = zeros(length(kappa_list), T, len_trial);
dists_X_GD = zeros(length(kappa_list), T, len_trial);

for trial = 1:len_trial
    freq_seed = ceil(n1*rand(r, 1))/n1;
    omega_seed = rand(nd, 1);
    w = zeros(nd, 1);
    for k = 1:nd
        w(k) = min([k, n1+n2-k, n1, n2]); %length of skew-diagonals
    end    
    for i_kappa = 1:length(kappa_list)
        kappa = kappa_list(i_kappa);
        sigma_star = linspace(kappa, 1, r);
        x_star = exp(2*pi*1i * (0:(nd-1))' * freq_seed') * sigma_star' / sqrt(n1*n2);
        omega = omega_seed < p;
        y = omega.*x_star;
        X_star = HankelLift(x_star, n1, n2);
        Y = HankelLift(y, n1, n2);
        [U0, Sigma0, V0] = svds(Y/p, r);
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
            x = HankelProj(L, R, nd, w);
            z = (omega.*x-y)/p - x;
            [Lz, Rz] = HankelMul(z, L, R, nd, n1, n2); 
            L_plus = L - eta*(L + Lz/(R'*R + eps('double')*eye(r)));
            R_plus = R - eta*(R + Rz/(L'*L + eps('double')*eye(r)));
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
            x = HankelProj(L, R, nd, w);
            z = (omega.*x-y)/p - x;
            [Lz, Rz] = HankelMul(z, L, R, nd, n1, n2);
            L_plus = L - eta/sigma_star(1)*(L*(R'*R) + Lz);
            R_plus = R - eta/sigma_star(1)*(R*(L'*L) + Rz);
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
    dists = dists(dists > 1e-12);
    T_subs = (2*i_kappa):10:length(dists);
    semilogy(T_subs, dists(T_subs), 'Color', clrs{1}, 'Marker', mks{i_kappa}, 'MarkerSize', 9);
    hold on; grid on;
    lgds{end+1} = sprintf('$\\mathrm{ScaledGD}~\\kappa=%d$', kappa);
end
for i_kappa = 1:length(kappa_list)
    kappa = kappa_list(i_kappa);
    dists = mean(dists_X_GD(i_kappa, :, :), 3);
    dists = dists(dists > 1e-12);
    T_subs = 10:10:length(dists);
    semilogy(T_subs, dists(T_subs), 'Color', clrs{2}, 'Marker', mks{i_kappa}, 'MarkerSize', 9);
    hold on; grid on;
    lgds{end+1} = sprintf('$\\mathrm{VanillaGD}~\\kappa=%d$', kappa);
end
xlabel('Iteration count');
ylabel('Relative error');
legend(lgds, 'Location', 'Best', 'Interpreter', 'latex');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
fig_name = sprintf('HankelMC_n=%d_r=%d_p=%g',n1,r,p);

function X = HankelLift(x, n1, n2)
% Lift the vector x as a Hankel matrix X = H[x]
    X = hankel(x);
    X = X(1:n1, 1:n2);  
end

function x = HankelProj(L, R, nd, w)
% Calculate the Hankel projection H[x] = H[L*R']
    x = sum(ifft(fft(L, nd, 1) .* fft(conj(R), nd, 1), nd, 1), 2) ./ w;
end

function [Lx, Rx] = HankelMul(x, L, R, nd, n1, n2)
% Calculate the Hankel multiplication Lx = H[x]*R, Rx = H[x]'*L
    Lx = ifft(bsxfun(@times, fft(R(n2:-1:1, :), nd, 1), fft(x)), nd, 1);
    Lx = Lx(n2:nd, :);  
    Rx = ifft(bsxfun(@times, fft(L(n1:-1:1, :), nd, 1), fft(conj(x))), nd, 1);
    Rx = Rx(n1:nd, :);
end