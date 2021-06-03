function M = real_matrix_sampling(A,eps)

M.exact = A;
[n_1,n_2] = size(M.exact);
m = round(eps*nnz(A));
perm = randperm(nnz(A));
gamma = find(A); %position of nonzero entries in A
M.ind_Omega = gamma(perm(1:m));
M.ind_Omega = sort(M.ind_Omega)';
[M.row,M.col] = ind2sub([n_1,n_2],M.ind_Omega);
M.values_Omega = M.exact(M.ind_Omega);
[M.size1,M.size2] = size(M.exact);
M.sparse = sparse(M.row,M.col,M.values_Omega,M.size1,M.size2);
M.proj  = @(X) X(M.ind_Omega)';
M.projT = @(y) sparse(M.row,M.col,y,M.size1,M.size2);