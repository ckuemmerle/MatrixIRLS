function [Mout, Out] = ScaledASD(m,n,r,Omega,data,start,opts)
%
% Implementation of matrix completion algorithm Scaled Alternating Steepest
% Descent ('ScaledASD') of [1]. 
% Code by Ke Wei.
% =========================================================================
% [1] Jared Tanner and Ke Wei. Low rank matrix completion by alternating 
% steepest descent methods. Applied and Computational Harmonic Analysis, 
% 40(2):417?429, 2016.
% =========================================================================
% Minor modifications by Christian Kuemmerle:
% - save intermediate iterates if opts.saveiterates == 1.
% - save timing information.

reltol = opts.rel_res_tol;
maxiter = opts.maxit;
verbosity = opts.verbosity;
rate_limit = 1-opts.rel_res_change_tol;
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    saveiterates = 1;
    Mout = cell(1,maxiter);
else
    saveiterates = 0;
end
relres = reltol * norm(data);

identity = eye(r);

p = length(data);

% creat a sparse matrix
[I, J] = ind2sub([m,n],Omega);
diff_on_omega_matrix = sparse(I,J,data,m,n,p);

if ~isempty(start)
  X = start.L;
  Y = start.R;
  if verbosity
    fprintf('initial point provided!\n')
  end
else
  X = randn(m,r);
  Y = randn(r,n);
  if verbosity
    fprintf('no initial point!\n')
  end
end

clear opts;
clear start;

Xt = X';
diff_on_omega = data'-partXY(Xt,Y,I,J,p);
res = norm(diff_on_omega);

iter = 1;

itres = zeros(maxiter,1);
itres(iter) = res;
time = zeros(1,maxiter);

conv_rate = 0;

tic
% iteration
while iter <= maxiter &&  res >=  relres && conv_rate <= rate_limit    
    % gradient for X
    updateSval(diff_on_omega_matrix,diff_on_omega,p);
    grad_X = diff_on_omega_matrix*Y'; % m*r
    
    % scaled gradient and stepsize for X
    scale = (Y*Y')\identity;
    dx = grad_X * scale; % m*r
    % dx = grad_X * diag(1./diag(Y*Y')); % scaled by diagonal entries
    delta_XY = partXY(dx',Y,I,J,p);
    tx = trace(dx'*grad_X)/norm(delta_XY)^2;
    
    % update X
    X = X + tx*dx;
    
    diff_on_omega = diff_on_omega-tx*delta_XY;
    
    % gradient for Y
    updateSval(diff_on_omega_matrix,diff_on_omega,p);
    Xt = X';
    grad_Y = Xt*diff_on_omega_matrix;
    
    % scaled gradient and stepsize for Y
    scale = (Xt*X)\identity;
    dy = scale*grad_Y;
    % dy = diag(1./diag(Xt*X))*grad_Y; % scaled by diagonal entries
    delta_XY = partXY(Xt,dy,I,J,p);
    ty = trace(dy*grad_Y')/norm(delta_XY)^2;
    
    % update Y
    Y = Y + ty*dy;

    diff_on_omega = diff_on_omega-ty*delta_XY;
    res = norm(diff_on_omega);
    
    % only for output
    if saveiterates
        Mout{iter}={X,Y'};
    end
    %
    
    iter = iter + 1;   
    itres(iter) = res;
    time(iter) = toc;
    conv_rate = (itres(iter)/itres(max(1,iter-15)))^(1/min(15,iter-1));
    
end
if saveiterates
    Mout = Mout(1:iter-1);
else
    Mout = {{X,Y'}};
end

Out.itrelres = itres(1:iter)/norm(data);
Out.iter = iter-1;
Out.reschg = abs(1-conv_rate);
Out.time = time(1:Out.iter);

end
