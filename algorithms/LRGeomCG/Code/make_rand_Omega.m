function [Omega, Omega_sparse] = make_rand_Omega(m,n,nb_samples)
% Make a random sampling using specified number of samples.
% Output is 
%    * Omega as a vector of (vectorized) indices.
%    * Omega_sparse as a sparse matrix (optional output).

if nb_samples<0 
  error('nb_samples has to be a fraction >=0')
end

if nb_samples<max(m,n)
  warning('You probably want to sample at least max(m,n) samples')
end

if nb_samples>m*n
  warning('Nb_samples is larger than max value (m*n). Setting to m*n although it will not result in an interesting matrix completion problem...')
  nb_samples = m*n;
end


Omega = randsampling(m,n,nb_samples);

if nargout==2
  [Omega_i, Omega_j] = ind2sub([m,n], Omega);
  Omega_sparse = sparse(Omega_i, Omega_j, 1, m, n, nb_samples);
end


end


function omega = randsampling(n_1,n_2,m)

omega = ceil(rand(m, 1) * n_1 * n_2);
omega = unique(omega);
while length(omega) < m    
    omega = [omega; ceil(rand(m-length(omega), 1)*n_1*n_2);];
    omega = unique(omega);
end
omega = omega(1:m);
omega = sort(omega);

end
