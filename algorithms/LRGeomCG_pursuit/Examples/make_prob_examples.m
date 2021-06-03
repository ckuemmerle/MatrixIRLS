function prob = make_prob(L,R,Omega,k, Gamma)
% L*R' is the exact data
% Omega is the sampling set; a sparse matrix.
% k is the rank to reconstruct
% Gamma (optional) is the testing set

if nargin==4
    Gamma = [];
end
prob.n1 = size(L,1); 
prob.n2 = size(R,1); 
prob.r = k;
if prob.r >= min(prob.n1, prob.n2)
    error('The rank of the matrices in the reconstruction manifold should be stricly smaller than number of rows and columns.')
end



prob.Omega = Omega;
[prob.Omega_i, prob.Omega_j] = ind2sub([prob.n1,prob.n2], Omega);
prob.m = length(Omega);

prob.data = partXY(L', R',prob.Omega_i, prob.Omega_j,prob.m)'; 
prob.temp_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data*1,prob.n1,prob.n2,prob.m);
  
% regularization parameter
%prob.mu = 1e-14;
prob.mu = 0;
prob.use_blas = 1;

if ~isempty(Gamma)    
    prob.has_testing = true;
    prob.Gamma = Gamma;
    [prob.Gamma_i, prob.Gamma_j] = ind2sub([prob.n1,prob.n2], Gamma);
    prob.test_data = partXY(L', R',prob.Gamma_i, prob.Gamma_j,prob.m); 
else
    prob.has_testing = false;
end
