function start = default_init(func_name,m,n,r,Omega,data)

[i,j] = ind2sub([m,n],Omega);
p = length(Omega);

M_omega = sparse(i, j, data, m, n, p);

[U,S,V] = svds(M_omega, r);

if strcmp(func_name,'RGrad') || strcmp(func_name,'RRCG') ||  strcmp(func_name,'RCG') 
  start.U = U;
  start.sigma = diag(S);
  start.V = V;
end

if strcmp(func_name,'ASD') || strcmp(func_name,'ScaledASD')
  start.L = U*sqrt(S);
  start.R = sqrt(S)*V';
end



