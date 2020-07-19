function opts = opts_adapted(maxit,rel_res_tol,rel_res_change_tol,...
    verbosity,saveiterates)

% maximun number of iteration
opts.maxit = maxit;%5000;

% Tolerance on the relative l_2 error 
% on the sampling set Omega (i.e., residual)
opts.rel_res_tol = rel_res_tol;%1e-5;

% Tolerance for detection of stagnation. 
opts.rel_res_change_tol = rel_res_change_tol;%1e-4;

% Verbosity 1 is chatty.
opts.verbosity = verbosity;%0;

% Saveiterates 1 means that all intermediate iterates will be saved.
opts.saveiterates = saveiterates;
