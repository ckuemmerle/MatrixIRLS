function opts = default_opts_pursuit()

opts.maxit = 100;

% Tolerance on the Riemannian gradient of the objective function
opts.abs_grad_tol = 0;
opts.rel_grad_tol = 0;
opts.abs_f_tol = 0;
opts.rel_f_tol = 0;


% Tolerance for detection of stagnation. 
opts.rel_tol_change_x = 0;
opts.rel_tol_change_res = 1e-4;

% Verbosity 2 is very chatty.
opts.verbosity = 1;

% Strong Wolfe is needed in theory but not usually in practice and is only a little
% slower
opts.strong_wolfe = true; 


opts.opt_S = 0; %very expensive

opts.rel_grad_decrease_factor = 1e-4;
opts.rand_pca = true;
opts.differentiated_transport = false;
