function x = make_start_x(prob, sparse_svd)
% Compute a truncated SVD for starting guess.
if nargin==1
    sparse_svd = true;
end

p = length(prob.Omega)/prob.n1/prob.n2;


M_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data,prob.n1,prob.n2,prob.m);
if sparse_svd
    opts.tol = 1e-5;
    [U,S,V] = svds(M_omega/p, prob.r, 'L', opts);
else
    [U,S,V] = svd(full(M_omega)/p);
    U = U(:,1:prob.r);
    S = S(1:prob.r,1:prob.r);
    V = V(:,1:prob.r);
end

x.sigma = [];
x.V = V;
x.S = S;
x.U = U;

x = prepx(prob, x);


