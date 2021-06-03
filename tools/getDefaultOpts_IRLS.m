function opts = getDefaultOpts_IRLS
% =========================================================================
% Reference:
% [1] C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples", ICML 2021.
% =========================================================================
opts.p              = 0;  % Schatten-p parameter (p close to 0 for very non-convex objective, p=1 for convex)
opts.N0             = 200;  % max. number of outer iterations (second order method like IRLS)
opts.N0_inner       = 200;  % max. number of inner iterations (second order method like IRLS)
opts.N0_firstorder  = 5000; % max. number of iteration for a first order method
opts.N_SVD          = 20;  % max. number of iterations for power method-type solver for partial SVD (such as bksvd)
opts.tol            = 1e-9;  % stopping criterion, lower bound on relative change of Frobenius norm
opts.tol_CG         = 1e-5; % tolerance for stopping criterion of CG (bound on relative residual norm)
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
opts.tangent_para   = 'extrinsic'; % decides what kind of computational representation of the tangent space is used: 'intrinsic' or 'extrinsic'
end

    