function small_example_decay()

randn('state',0); rand('state',0);

m = 1000; n = 1000; % dimensions
k = 40; % rank
% Relative oversampling factor 
% OS=1 is minimum, 2 is difficult, 3 is OKish, 4+ is easy.
OS = 5;


% random factors
L0 = randn(m, k); [L0,~] = qr(L0,0);
R = randn(n, k); [R,~] = qr(R,0);
S = diag(logspace(0,-15,k)); L = L0*S;
dof = k*(m+n-k);
rank_reconstruct = 10;

% make random sampling, problem and initial guess
samples = floor(OS * dof);
Omega = make_rand_Omega(m,n,samples);
prob = make_prob(L,R,Omega,rank_reconstruct); % <- you can choose another rank here

options = default_opts();
options.maxit = 100;
x0 = make_start_x(prob);

t=tic;
[Xcg,hist] = LRGeomCG(prob,options,x0);
out_time = toc(t)
Xcg = Xcg.U * diag(Xcg.sigma) * Xcg.V';

semilogy(hist(:,1),'rx')
hold on
semilogy(hist(:,2),'bx')
legend('Rel grad', 'Rel error on Omega')

X = L*R';
error_cg = norm(Xcg-X,'fro')/norm(X,'fro')
Xbest = L0(:,1:rank_reconstruct)*S(1:rank_reconstruct,1:rank_reconstruct)*R(:,1:rank_reconstruct)';
error_bestrank = norm(Xbest-X,'fro')/norm(X,'fro')

diff_cg_bestrank = norm(Xbest-Xcg,'fro')/norm(Xcg,'fro')


