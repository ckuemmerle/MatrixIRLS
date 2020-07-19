% Test code for the RTRMC/RCGMC algorithm (low-rank matrix completion).
% Algorithm by Nicolas Boumal and P.-A. Absil.
% Code by Nicolas Boumal, UCLouvain, April 2, 2013.
%
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
% http://sites.uclouvain.be/absil/RTRMC/
%
% See also: rtrmc buildproblem initialguess

clear all;
close all;
clc;

if exist('manoptsolve', 'file') ~= 2
    cd 'Manopt_2.0';
    importmanopt;
    cd ..;
    warning('manopt:version', ['Using Manopt 2.0. This might not be ' ...
             'the most up-to-date version of Manopt. Please consider going to ' ...
             'http://www.manopt.org to obtain the latest version.\n\n' ...
             'Manopt 2.0 was added to your Matlab path.\n\n']);
end

% If this script fails, try executing installrtrmc.m, to compile the mex files.
% This will not be necessary if launching the present script works out of the box,
% which there is a decent chance it should.

%% Problem instance generation

% Dimensions of the test problem
m = 5000;                              % number of rows
n = 5000;                              % number of columns
r = 5;                                 % rank
k = floor(2*r*(m+n-r));                       % number of known entries

% Generate an m-by-n matrix of rank true_rank in factored form: A*B
true_rank = r;
A = randn(m, true_rank)/true_rank.^.25;
B = randn(true_rank, n)/true_rank.^.25;

% Pick k (or about k) entries uniformly at random
[I, J, k] = randmask(m, n, k);

% Compute the values of AB at these entries
% (this is a C-Mex function)
X = spmaskmult(A, B, I, J);

% Define the confidence we have in each measurement X(i)
% IMPORTANT:
% If you wish to use non-uniform weights (so, C is not a constant vector),
% you need the latest compiled version of spbuildmatrix.c. You can obtain
% it by executing:
%   mex -lmwlapack -lmwblas -largeArrayDims spbuildmatrix.c
% You do not need to do that if you are using .mexw64 or .mexmaci64 files
% (they are already up to date).
% If you compile the file for a new platform, it would be great if you sent
% us the output of that compilation, for inclusion in the package:
% nicolasboumal@gmail.com.
% Thanks!
C = ones(size(X));

% Add noise if desired
noisestd = 0;
X = X + noisestd*randn(size(X));


%% Feeding the problem instance to RTRMC

% Pick a value for lambda, the regularization parameter
lambda = 0;


perm = randperm(k);
I = I(perm);
J = J(perm);
X = X(perm);
C = C(perm);


% Build a problem structure
problem = buildproblem(I, J, X, C, m, n, r, lambda);


% Compute an initial guess
initstart = tic;
U0 = initialguess(problem);
inittime = toc(initstart);

% [Optional] If we want to track the evolution of the RMSE as RTRMC
% iterates, we can do so by specifying the exact solution in factored form
% in the problem structure and asking RTRMC to compute the RMSE in the
% options structure. See the subfunction computeRMSE in rtrmc.m to see how
% to compute the RMSE on a test set if the whole matrix is not available in
% a factorized A*B form (which is typical in actual applications).
problem.A = A;
problem.B = B;

% Setup the options for the RTRMC algorithm
% These are the algorithms shown in the papers:
%  RTRMC 2  : method = 'rtr', order = 2, precon = false
%  RTRMC 2p : method = 'rtr', order = 2, precon = true
%  RTRMC 1  : method = 'rtr', order = 1, precon = false
%  RCGMC    : method = 'cg',  order = 1, precon = false
%  RCGMCp   : method = 'cg',  order = 1, precon = true
opts.method = 'rtr';     % 'cg' or 'rtr', to choose the optimization algorithm
opts.order = 2;          % for rtr only: 2 if Hessian can be used, 1 otherwise
opts.precon = true;      % with or without preconditioner
opts.maxiter = 1000;      % stopping criterion on the number of iterations
opts.maxinner = 50;      % for rtr only : maximum number of inner iterations
opts.tolgradnorm = 1e-8; % stopping criterion on the norm of the gradient
opts.verbosity = 2;      % how much information to display during iterations
opts.computeRMSE = true; % set to true if RMSE is to be computed at each step


opts.saveiterates = 1;

% Call the algorithm here. The outputs are U and W such that U is
% orthonormal and the product U*W is the matrix estimation. The third
% output, stats, is a structure array containing lots of information about
% the iterations made by the optimization algorithm. In particular, timing
% information should be retrieved from stats, as in the code below.
[U, W, stats] = rtrmc(problem, opts, U0);


time = inittime + [stats.time];
rmse = [stats.RMSE];

switch opts.method
    case 'rtr'
        semilogy(time, rmse, '.-');
    case 'cg'
        semilogy(time, rmse, '.-');
        S = [stats.linesearch];
        ls_evals = [S.costevals];
        for i = 1 : length(ls_evals)
            text(time(i), rmse(i)*1.20, num2str(ls_evals(i)));
        end
        title(sprintf('Numbers indicate line search cost evaluations.'));
end
xlabel('Time [s]');
ylabel('RMSE');
