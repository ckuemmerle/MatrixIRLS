function small_example()

randn('state',0); rand('state',0);

m = 1000; n = 1000; % dimensions
k = 40; % rank
% Relative oversampling factor 
% OS=1 is minimum, 2 is difficult, 3 is OKish, 4+ is easy.
OS = 2;


% random factors
L = randn(m, k); 
R = randn(n, k); 
dof = k*(m+n-k);

% make random sampling, problem and initial guess
samples = floor(OS * dof);
Omega = make_rand_Omega(m,n,samples);
prob = make_prob(L,R,Omega,k); % <- you can choose another rank here

options = default_opts();
x0 = make_start_x(prob);

t=tic;
[Xcg,hist] = LRGeomCG(prob,options,x0);
out_time = toc(t)

semilogy(hist(:,1),'rx')
hold on
semilogy(hist(:,2),'bx')
legend('Rel grad', 'Rel error on Omega')

