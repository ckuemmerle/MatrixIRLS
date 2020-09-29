function opts = default_opts_R3MC
% Sets the default parameters for the low-rank matrix completion method
% 'R3MC' based on "Test_R3MC.m" by Bamdev Mishra, 2013.
% =========================================================================
% Author: Christian Kuemmerle, Johns Hopkins University, kuemmerle@jhu.edu,
% 2020.
% =========================================================================

tol = 1e-10; % Tolerance; Look at the R3MC file to stop R3MC in a different way
verbose = true; % Show output
N0 = 500; % Maximum number of iterations allowed
ls_maxiter = 50; % Maximum linesearch allowed at each iteration
beta_type = 'F-R'; % Algorithm type. Other choices include H-S, F-R and off, P-R
linearized_linesearch = true; % Whether to use linearized stepsize search
accel_linesearch = true; % Whether to use accelerated linesearch
compute_predictions = false; % If interested in compute recovery at each iteration
random_initialization = false; % How to initialize the iterates
ORTH_VALUE = Inf; % Inf gives the best results for R3MC
sigma_armijo = 1e-4; % for Armijo linesearch. Keep a low valaue while going for linearized stpesize search
if ~linearized_linesearch
    sigma_armijo = 0.5;
end

opts.tol = tol;
opts.N0 = N0;
opts.sigma_armijo = sigma_armijo;
opts.compute_predictions = compute_predictions;
opts.beta_type = beta_type;
opts.linearized_linesearch =  linearized_linesearch;
opts.verbose = verbose;
opts.orth_value = ORTH_VALUE;
opts.accel_linesearch = accel_linesearch;
opts.random_initialization = random_initialization;
opts.accel_linesearch = true;
opts.ls_maxiter = ls_maxiter;
opts.gamma = 0;
opts.rank_updates = false;
opts.start_rank = 1;
end

