function U0 = initialguess(problem)
% FUNCTION U0 = INITIALGUESS(PROBLEM)
%
% Generates an initial guess of the column space of X, the matrix to be
% recovered, based on the observed entries of X and the mask pattern.
%
% Input:
%
% PROBLEM: A structure describing the low-rank matrix completion problem
%          to solve. Such a structure may be built using BUILDPROBLEM.
%
% Output:
%
% U0: An m-by-r orthonormal matrix spanning a column space that should be
%     closer to the true column space of the matrix to be recovered than a
%     simple random guess. U0 can be fed to RTRMC as initial guess to be
%     improved.
%
% Nicolas Boumal, UCLouvain, Sept. 3, 2012.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: buildproblem rtrmc

    m = problem.m;
    n = problem.n;
    k = problem.k;
    r = problem.r;
    I = problem.I;
    J = problem.J;
    X = problem.X;

    [U0, ~, ~] = svds(sparse(double(I), double(J), X, m, n, k), r);
    
end
