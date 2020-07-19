function opts = default_opts()

% maximun number of iteration
opts.maxit = 5000;

% Tolerance on the relative l_2 error 
% on the sampling set Omega (i.e., residual)
opts.rel_res_tol = 1e-5;

% Tolerance for detection of stagnation. 
opts.rel_res_change_tol = 1e-4;

% Verbosity 1 is chatty.
opts.verbosity = 0;
