function problem = buildproblem(I, J, X, C, m, n, r, lambda)
% PROBLEM = BUILDPROBLEM(I, J, X, C, M, N, R, LAMBDA)
%
% Build a structure describing a low-rank matrix completion problem, to be
% used with RTRMC.
%
% Inputs:
%
% I, J, X, C: Column vectors of equal length such that the k-th known entry
%             of the matrix to recover is at position (I(k), J(k)) and has
%             value X(k), with confidence C(k). If C is set to the empty
%             matrix, it will be set to ones(size(X)).
% 
% M, N, R:    The matrix to recover is M-by-N and assumed to be of rank R.
%             RTRMC works best if M <= N. If this is not the case, one
%             should build the problem as follows:
%             BUILDPROBLEM(J, I, X, C, N, M, R, LAMBDA),
%             i.e., build a problem structure for the transpose of the
%             matrix to recover, call RTRMC to solve the problem, which
%             yields factors U and W, then obtain the reconstructed matrix
%             as W'*U' instead of U*W.
%
% LAMBDA:     A small smoothing factor, set in accordance with the noise
%             level.
% IMPORTANT!  If LAMBDA is to be changed later on, one should
%             call BUILDPROBLEM to obtain a new structure and not simply
%             change the value of LAMBDA in the existing problem structure,
%             seen as LAMBDA is used to precompute a number of other fields
%             in the problem structure.
% 
% 
% Output:
%
% PROBLEM:    A structure to be passed to RTRMC to solve the low-rank
%             matrix completion problem.
%
% Nicolas Boumal, UCLouvain, Sept. 6, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% Modified on March 6, 2014 to prevent a former bug linked to the ordering
% of vectors I, J, X and C. The function will now return a proper problem
% structure for use with the RTRMC algorithms. Fields I, J, X and C in this
% problem structure may have been re-ordered (consistently).
%
% SEE ALSO: rtrmc initialguess

    if m > n
        warning('RTRMC:baddimensions', ...
         ['RTRMC is optimized to deal with matrices such that m <= n. ' ...
          'Please transpose your matrix.']);
    end
    
    if ~(strcmpi(class(I), 'uint32') && strcmpi(class(J), 'uint32'))
        I = uint32(I);
        J = uint32(J);
        % warning('RTRMC:uint32', ...
        %           'The index vectors I and J were converted to uint32.');
    end
    
    if isempty(C)
        C = ones(size(X));
    end
    
    %% Some verification of the consistency of the data
    assert(length(J) == length(I), 'I and J must have the same length.');
    assert(length(X) == length(I), 'I and X must have the same length.');
    assert(length(C) == length(I), 'I and C must have the same length.');
    assert(all(I >= 1) && all(I <= m), 'm is the number of rows, indexed by I.');
    assert(all(J >= 1) && all(J <= n), 'n is the number of columns, indexed by J.');
    assert(r <= min(m, n), 'The desired rank r may not be larger than min(m, n).');
    
    %% Reorder the vectors I, J, X and C such that their entries are listed
    %% in the same order as they would be if they were returned by calling
    %% the function find() on a sparse matrix containing these data.
    %% This is necessary because of the "hack" used in the C-Mex function
    %% setsparseentries to speed up Matlab processing of sparse matrices.
    %%  Added March 6, 2014.
    [~, order] = sort(sub2ind([m, n], I, J));
    I = I(order);
    J = J(order);
    X = X(order);
    C = C(order);

    %% Input argument data
    problem.I = I;
    problem.J = J;
    problem.X = X;
    problem.C = C;
    problem.m = m;
    problem.n = n;
    problem.r = r;
    problem.lambda = lambda;
    
    %% Dependent data we 'precompute'
    problem.k = length(I);
    problem.sr = problem.k/(m*n);
    problem.Chat = C.^2 - lambda.^2;
    
    
    % mask is a sparse matrix with the same sparsity structure as X, C,
    % Chat etc. We use it as a place holder to accelerate some computations
    % and memory manipulations that Matlab has a hard time with.
    problem.mask = sparse(double(problem.I), double(problem.J), ...
                          ones(problem.k, 1), m, n, problem.k);

    
    %% Identification of the problem structure
    problem.id = cputime;
    
    % This is a limitation due to the implementation of buildAchol.c,
    % it is not a mathematical limitation. This being said, it does not in
    % general make sense to have negative Chat entries, as this means that
    % some observed entries are less certain than the unobserved entries.
    % To make this right, simply make sure all entries in C are >= lambda.
    assert(all(problem.Chat >= 0), ...
        ['All entries in C (the weights) must be larger than lambda ' ...
         '(the regularization).']);

end
