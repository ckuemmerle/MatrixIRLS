function start = make_start_x(func_name,m,n,r,Omega,data)

[i,j] = ind2sub([m,n],Omega);
p = length(Omega);

M_omega = sparse(i, j, data, m, n, p);

[U,S,V] = svds(M_omega, r);

if strcmp(func_name,'CGIHT_Matrix') || strcmp(func_name,'NIHT_Matrix') 
  start.U = U;
  start.sigma = diag(S);
  start.V = V;
end

if strcmp(func_name,'ASD') || strcmp(func_name,'ScaledASD')
  start.L = U*S;
  start.R = V';
end



