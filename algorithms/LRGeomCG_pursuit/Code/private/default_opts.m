function opts = default_opts()

opts.maxit = 5000;

% Tolerance on the Riemannian gradient of the objective function
opts.abs_grad_tol = 0;
opts.rel_grad_tol = 1e-8;

% Tolerance on the l_2 error on the sampling set Omega
opts.abs_f_tol = 0;
opts.rel_f_tol = 1e-12;

% Tolerance for detection of stagnation. 
opts.rel_tol_change_x = 1e-12;
opts.rel_tol_change_res = 1e-4;

% Verbosity 2 is very chatty.
opts.verbosity = 1;

% Strong Wolfe is needed in theory but not in practice and is a little
% slower
opts.strong_wolfe = false; 

opts.rel_grad_decrease_factor = 0;

opts.differentiated_transport = false; % make no difference and is slower, but needed for convergence proof

