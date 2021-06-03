function [prob,x] = make_prob(n, rank, m, state)

if nargin==4
    randn('state', state);
    rand('state', state);
end

prob.n1 = n(1); prob.n2 = n(2); prob.r = rank;


prob.M_exact_A = randn(prob.n1,prob.r);
prob.M_exact_B = randn(prob.n2,prob.r);

%[prob.M_exact_A,Ra] = qr(prob.M_exact_A,0);
%[prob.M_exact_B,Rb] = qr(prob.M_exact_B,0);
%s = 10.^(-(1:rank));
%prob.M_exact_A = prob.M_exact_A*diag(s);

df = prob.r*(prob.n1+prob.n2-prob.r);
prob.m = m;

prob.Omega = randsample(prob.n1*prob.n2,prob.m);
prob.Omega = sort(prob.Omega);
%[temp,prob.Omega_indx] = sort(prob.Omega);
[prob.Omega_i, prob.Omega_j] = ind2sub([prob.n1,prob.n2], ...
    prob.Omega);


prob.data =  XonOmega(prob.M_exact_A, prob.M_exact_B, prob.Omega);
prob.M_Omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data,prob.n1,prob.n2,prob.m);
%full(prob.M_Omega(1:5,1:5))
%updateSparse(prob.M_Omega, prob.data, prob.Omega_indx);
%full(prob.M_Omega(1:5,1:5))
prob.temp_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data*1,prob.n1,prob.n2,prob.m);

m_df = prob.m / df;
m_n_n = prob.m / prob.n1 / prob.n2;

% check 
%T = zeros(n1,n2);
%M = prob.M_exact_A*prob.M_exact_B';
%T(prob.Omega) = M(prob.Omega);
%norm(prob.M_Omega - T,'fro')

k = rank;
x.U = randn(prob.n1,k); [x.U,temp] = qr(x.U,0);
x.sigma = sort(abs(randn(k,1)),'descend');
x.V = randn(prob.n2,k); [x.V,temp] = qr(x.V,0);
x = prepx(prob, x);
