function [Mout, Out] = CGIHT_Matrix(m,n,r,Omega,data,start,opts)
%
% Implementation of matrix completion algorithm Conjugate Gradient
% Iterative Hard Thresholding ('CGIHT') of [1]. 
% Code by Ke Wei.
% =========================================================================
% [1] Jeffrey Blanchard, Jared Tanner, and Ke Wei. CGIHT: Conjugate gradient 
% iterative hard thresholding for compressed sensing and matrix completion. 
% Information and Inference, 4(4):289?327, 2015.
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

% note that sparse(full) is necessary
% or else d will point to the same area as res_on_omega_matrix
d = sparse(full(res_on_omega_matrix));
Ad_proj = partXY(U',U'*d, I, J, p);

% iteration
while  ((res >= relres) && (iter <= maxiter)) % && (conv_rate < rate_limit)) 
    % compute alpha
    % <UU'*r, UU'*d> = <r, UU'*d> = <r, A.*(UU'*d)>
    % since res_on_omega is sparse
    % alpha = d_proj*res_on_omega/norm(d_proj)^2;

    alpha = Ad_proj*res_on_omega;
    tmp = norm(Ad_proj)^2;
    if abs(alpha) < 1000*tmp
      alpha = alpha/tmp;
    else
      alpha = 1;
    end
    
    % compute M+alpha*p_new
    M = M + alpha * d;
    
    [U,S,V] = svds(M,r);
    M = U*S*V'; % can be further accelerated here
    if saveiterates
        Mout{iter}={U*S,V};
    end
    
    res_on_omega = data - M(Omega); % p*1 vector
    updateSval(res_on_omega_matrix,res_on_omega,p);    
    res = norm(res_on_omega);

    % res_proj = A.*(UU'*res_on_omega) = A.*(U(U'*res_on_omega))
    % d_proj = A.*(UU'*d) = A.*(U(U'*d))
    Ares_proj = partXY(U', U'*res_on_omega_matrix, I, J, p);
    Ad_proj = partXY(U', U'*d, I, J, p);
    
    %beta = res_proj*d_proj'/norm(d_proj)^2; 
    beta = Ares_proj*Ad_proj';
    tmp = norm(Ad_proj)^2;

    if abs(beta) < 1000 * tmp
      beta = beta/tmp;
    else
      beta = 0;
    end
   
    d = res_on_omega_matrix - beta*d; 

    % update d_proj
    Ad_proj = Ares_proj - beta * Ad_proj;

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

