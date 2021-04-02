% Modification of LRGeomC
rng(9)

% Choose your example A,B,C below by uncommenting/commenting the correct
% lines.

options = default_opts_pursuit();
options.increase_with_residual = true;
options.strong_wolfe = false;

% % ========================================
% % A) Example of paper: Random matrix with exponential decaying sing values
% m = 1000; n = 1000; k = 26; % dimensions and rank
% OS = 2.5; % oversampling factor  OS=1 is minimum, 2 is difficult, 3 is OKish, 4+ is easy.
% prob = make_prob_randomLR(m,n,logspace(0,-10,k), k, OS, true);
% 
% options.rank_increase = 1; % choose 1 or 2
% options.rel_grad_decrease_factor = 1e-5; 
% options.max_rank = k+2;

% % ========================================
% % B) Example of paper: Bivariate function
% m = 1000; n = 1000; k = 20; % dimensions and rank
% OS = 2.5; % oversampling factor  OS=1 is minimum, 2 is difficult, 3 is OKish, 4+ is easy.
% fun_f = @(x,y)(1./(1+(x-y).^2)); prob = make_prob_function(fun_f, [0 1 0 1], m,n, k,OS,true);
% 
% options.rank_increase = 1; % choose 1 or 2
% options.rel_grad_decrease_factor = 1e-4; 
% options.max_rank = 18;

% ========================================
% C) New example of paper: Plateau of singular values. Will work less well but it is still surprisingly OK.
m = 1000; n = 1000; k = 26; % dimensions and rank
OS = 1.5; % oversampling factor  OS=1 is minimum, 2 is difficult, 3 is OKish, 4+ is easy.
prob = make_prob_randomLR(m,n,[ones(1,7) logspace(-1,-5,8) 1e-6*ones(1,k-7-8)], k, OS, true);

options.rank_increase = 2; % choose 1 or 2
options.rel_grad_decrease_factor = 1e-5; 
options.max_rank = k+2;
options.saveiterates = 1;

[x,hist,stats,Xr,outs] = LRGeomCG_pursuit(prob, options);

%% Plot
figure(1)
semilogy(hist(:,1),'ro')
hold on
semilogy(hist(:,2),'bo')
semilogy(hist(:,5),'gd')
legend('Rel grad', 'Rel training error', 'Rel testing error')
% put vertical lines
for ii=1:length(stats.chg_ranks)
    semilogy([stats.chg_ranks(ii) stats.chg_ranks(ii)], [min(hist(:,1)) 1], 'k--', 'HandleVisibility','off')
end
ylim([1e-15 1])
xlabel('Inner iterations')
