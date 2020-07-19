function opts = getDefaultOpts_IRLS

opts.p              = 0;  % Schatten-p parameter (p close to 0 for very non-convex objective, p=1 for convex)
opts.N0             = 200;  % max. number of outer iterations (second order method like IRLS)
opts.N0_inner       = 50;  % max. number of inner iterations (second order method like IRLS)
opts.N0_firstorder  = 5000; % max. number of iteration for a first order method
opts.tol            = 1e-9;  % stopping criterion, lower bound on relative change of Frobenius norm
opts.tol_CG_fac     = 1e-5; % multiplicative factor such that tol_CG = eps_c * tol_CG_fac.
opts.epsmin         = 1e-15;  % minimal value for epsilon smoothing
opts.use_min        = 1;
opts.type_mean      = 'optimal';
opts.mode_linsolve  = 'tangspace';
opts.qmean_para     = min(opts.p/(opts.p-2),0);
opts.increase_antisymmetricweights = 0;
opts.mode_eps       = 'oracle_model_order';%'classical_r' % iter_diff'
opts.objective      = 'objective_thesis';
opts.tracking = 0;
opts.adaptive_cg    = 0;
opts.verbose        = 1;
end

    