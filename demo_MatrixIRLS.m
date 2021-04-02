% This is a minimal example demonstrating how to use the implementation of 
% Matrix Iteratively Reweighted Least Squares ('MatrixIRLS'), cf.
% [1,2], to complete low-rank matrices.
% =========================================================================
% References:
% [1] C. Kuemmerle, C. M. Verdun, "Escaping Saddle Points in 
% Ill-Conditioned Matrix Completion with a Scalable Second Order Method", 
% ICML 2020 Workshop "Beyond First Order Methods in ML Systems".
%
% [2] C. Kuemmerle, "Understanding and Enhancing Data Recovery Algorithms From
% Noise-Blind Sparse Recovery to Reweighted Methods for Low-Rank Matrix 
% Optimization", Ph.D. dissertation, Technical University of Munich, 2019,
% Chapter 2.
% Available at: https://mediatum.ub.tum.de/doc/1521436/1521436.pdf
% =========================================================================
% Author: Christian Kuemmerle, 2020.
%% Set parameters
% Number of rows of the matrix that defines the problem:
rng('shuffle');
d1 = 1000; 
% Number of columns:
d2 = 1000; 
% Rank:
r = 10;
% Number of degrees of freedom of the setting:
df_LR = @(rr) rr*(d1 + d2 - rr);
df_LR_val=df_LR(r);
% Oversampling factor: % Choose between 1-1.5 for very hard problems
% (sometimes with no unique solution), or > 2 for easier problems.
oversampling = 2;
% Number of measurements:
m = floor(min(oversampling*df_LR_val,d1*d2));
%% Sample the measurement matrix Phi (revealed entries correspond to non-zero entries) 
[Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,50);
%% Sample the ground truth matrix X0 and measured entries y
modeX0      = 2;
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,0);
X0 = {U0,V0}; % factorized representation of matrix X0= U0*V0' to be completed
X0_full = U0*V0';
y = X0_full(Omega); % vector with observed entries

%% Choose algorithm, set algorithmic options
opts = getDefaultOpts_IRLS;
opts.p = 0; % (non-)convexity parameters p for IRLS
% p = 0: sum of log objective
% 0 < p < 1: Schatten-p quasi-norm.
% p = 1: Nuclear norm.
opts.N0_inner = 200;
opts.saveiterates = 1;
opts.verbose = 1;
opts.tol = 1e-12;
opts.tangent_para = 'extrinsic';
%% Run MatrixIRLS
prob.d1     = d1;
prob.d2     = d2;
prob.r      = r;
prob.Phi    = Phi;
prob.y      = y;
lambda = 0;
[X_c,outs] = MatrixIRLS(prob,lambda,opts);

%% Get reconstructed matrix
X_c_full = get_densemat_from_compact(X_c,Phi);
%% Calculate the error of iterates to ground truth X0
alg_names = {'MatrixIRLS'};
Xr = cell(1,1);
if opts.saveiterates
    for k=1:outs.N 
        Xr{1}{k} = outs.X{k};
    end
else
    Xr{1}{1} = X_c;
end
[error_fro_rel,error_fro] = get_frob_errors(Xr,X0,Phi,alg_names,...
    'full',1);
%% Visualization of Frobenius errors of iterates
visualize_errorcurves_combined(error_fro_rel,alg_names);
plot_times_errors({outs},error_fro_rel,alg_names);