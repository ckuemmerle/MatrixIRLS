function opts = default_opts_Riem(N0,tol)

opts.maxit = N0;

% Tolerance on the Riemannian gradient of the objective function
opts.abs_grad_tol = min(1e-4,tol);%1e-7;
opts.rel_grad_tol = tol;%1e-10;

% Tolerance on the l_2 error on the sampling set Omega
opts.abs_f_tol = min(1e-4,tol);%1e-7;
opts.rel_f_tol = tol;%1e-10;

% Tolerance for detection of stagnation. 
opts.stagnation_detection = 1;  % <--- cheap test but does take some time
opts.rel_tol_change_x = tol;%1e-7;
opts.rel_tol_change_res = tol;%1e-6;

% Verbosity 2 is very chatty.
opts.verbosity = 1;

