function[U, Sigma, V] = svd_rank_1_update(U,Sigma,V,u,v)
%
% SVD of USigmaV' + u*v'
% Refer to paper
% Fast Low-Rank Modifications of the Thin Singular Value Decomposition, Matthew Brand, 2006.

r = size(U,2);
m = U'*u;
p = u - U*m;
R_u = sqrt(p'*p);
P = p/R_u;

n = V'*v;
q = v - V*n;
R_v = sqrt(q'*q);
Q = q/R_v;

% Update
U(:,r+1) = P;
V(:,r+1) = Q;
Sigma(:,r+1) = 0;
Sigma(r+1,:) = 0;

L = Sigma + [m ; R_u]*[n ; R_v]';

[A S B] = svd(L);

U = U*A;
V = V*B;
Sigma = S;
