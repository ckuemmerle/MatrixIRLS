function opts = default_opts(maxit,rel_res_tol,rel_res_chg_tol)
% maximun number of iteration
opts.maxit = maxit;

% tolerance on the relative residual
opts.rel_res_tol = rel_res_tol;

% tolerance for detection of stagnation. 
opts.rel_res_chg_tol = rel_res_chg_tol;

