function [Mout, Out] = NIHT_Matrix(m,n,r,Omega,data,start,opts)
%
% Implementation of matrix completion algorithm Normalized Iterative Hard 
% Thresholding ('NIHT') of [1]. 
% Code by Ke Wei.
% =========================================================================
% [1] Jared Tanner and Ke Wei. Normalized iterative hard thresholding 
% for matrix completion. SIAM J. Scientfic Comput., 35, S104?S125, 2013.
% =========================================================================
% Minor modifications by Christian Kuemmerle:
% - save intermediate iterates if opts.saveiterates == 1.
% - save timing information.

reltol = opts.rel_res_tol; 
maxiter = opts.maxit;
rate_limit = 1-opts.rel_res_change_tol;
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    saveiterates = 1;
    Mout=cell(1,maxiter);
else
    saveiterates = 0;
end
relres = reltol*norm(data);
itres = zeros(maxiter,1);
conv_rate = 0;

U = start.U;
S = diag(start.sigma);
V = start.V;

clear opts;
clear start;

p = length(data);
[I,J] = ind2sub([m,n],Omega);
res_on_omega_matrix = sparse(I,J,data,m,n,p);

M = U*S*V';
res_on_omega = data - M(Omega); % p*1 vector
updateSval(res_on_omega_matrix,res_on_omega,p);
res = norm(res_on_omega);

iter = 1;
itres(iter) = res;
time = zeros(1, maxiter);
tic

% iteration
while  ((res >= relres) && (iter <= maxiter))
    % compute alpha
    % <UU'*r, UU'*d> = <r, UU'*d> = <r, A.*(UU'*d)>
    % since res_on_omega is sparse
    % alpha = Ares_proj*res_on_omega/norm(d_proj)^2;
    Ares_proj = partXY(U', U'*res_on_omega_matrix, I, J, p);
    alpha = Ares_proj*res_on_omega;
    tmp = norm(Ares_proj)^2;
    if abs(alpha) < 1000*tmp
      alpha = alpha/tmp;
    else
      alpha = 1;
    end
    
    % compute M+alpha*p_new
    M = M + alpha * res_on_omega_matrix;
    
    [U,S,V] = svds(M,r);
    M = U*S*V'; % can be further accelerated here
    if saveiterates
        Mout{iter}={U*S,V};
    end

    res_on_omega = data - M(Omega); % p*1 vector
    updateSval(res_on_omega_matrix,res_on_omega,p);    
    res = norm(res_on_omega);

    iter = iter + 1;
    
    itres(iter) = res;
    conv_rate = (itres(iter)/itres(max(1,iter-15)))^(1/min(15,iter-1));
    time(iter) = toc;
    if iter >= 500 && conv_rate >= rate_limit
        break;
    end
end
iter=iter-1;
if saveiterates
    Mout = Mout(1:iter);
else
    Mout = {{U*S,V}};
end

Out.itrelres = itres(1:iter)/norm(data);
Out.iter = iter;
Out.reschg = abs(1-conv_rate);
Out.time = time(1:iter);

