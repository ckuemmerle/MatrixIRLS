function [varargout] = eigens(A,k,its,l)
%EIGENS  Eigendecomposition of a SELF-ADJOINT matrix
%
%
%   [V,D] = EIGENS(A)  constructs a nearly optimal rank-6 approximation
%           VDV' to A, using 4 normalized power iterations, with block
%           size 6+2=8, started with an n x 8 random matrix, when A is
%           n x n; the reference below explains "nearly optimal." When
%           A is the only input to EIGENS, the dimension n of A must be
%           >= 6.
%
%   [V,D] = EIGENS(A,k)  constructs a nearly optimal rank-k approx.
%           VDV' to A, using 4 normalized power iterations, with block
%           size k+2, started with an n x (k+2) random matrix, when A
%           is n x n; the reference below explains "nearly optimal."
%           k must be a positive integer <= the dimension n of A.
%
%   [V,D] = EIGENS(A,k,its)  constructs a nearly optimal rank-k approx.
%           VDV' to A, using its normalized power iterations, with block
%           size k+2, started with an n x (k+2) random matrix, when A
%           is n x n; the reference below explains "nearly optimal."
%           k must be a positive integer <= the dimension n of A, and
%           its must be a nonnegative integer.
%
%   [V,D] = EIGENS(A,k,its,l)  constructs a nearly optimal rank-k
%           approximation VDV' to A, using its normalized power
%           iterations, with block size l, started with an n x l random
%           matrix, when A is n x n; the reference below explains
%           "nearly optimal." k must be a positive integer <= the
%           dimension n of A, its must be a nonnegative integer, and l
%           must be a positive integer >= k.
%
%
%   Replacing "[V,D]" with "lambda" in any of the above produces the
%   vector lambda = diag(D) rather than both V and D.
%
%
%   The low-rank approximation VDV' comes in the form of an
%   eigendecomposition -- the columns of V are orthonormal and D is a
%   real diagonal matrix such that the absolute values of its diagonal
%   entries are nonincreasing.
%   V is n x k and D is k x k, when A is n x n.
%
%   Increasing its or l improves the accuracy of the approximation VDV';
%   the reference below describes how the accuracy depends on its and l.
%
%
%   Note: THE MATRIX A MUST BE SELF-ADJOINT.
%
%   Note: EIGENS invokes RAND. To obtain repeatable results, reset the
%         seed for the pseudorandom number generator.
%
%   Note: The user may ascertain the accuracy of the approximation VDV'
%         to A by invoking DIFFSNORMS(A,V,D).
%
%
%   inputs (the first is required):
%   A -- matrix being approximated
%   k -- rank of the approximation being constructed;
%        k must be a positive integer <= the dimension of A, and
%        defaults to 6
%   its -- number of normalized power iterations to conduct;
%          its must be a nonnegative integer, and defaults to 4
%   l -- block size of the normalized power iterations;
%        l must be a positive integer >= k, and defaults to k+2
%
%   outputs (produces either the first two, V and D, or just lambda):
%   V -- n x k matrix in the rank-k approximation VDV' to A, where A is
%        n x n
%   D -- k x k diagonal matrix in the rank-k approximation VDV' to A,
%        such that the diagonal entries are real and their absolute
%        values are nonincreasing
%   lambda -- k x 1 vector equal to diag(D)
%
%
%   Example:
%     A = rand(1000,2);
%     A = A*A';
%     A = A/normest(A);
%     [V,D] = eigens(A,2);
%     diffsnorms(A,V,D)
%
%     This example produces a rank-2 approximation VDV' to A such that
%     the columns of V are orthonormal and D is a diagonal matrix whose
%     diagonal entries are real and their absolute values are
%     nonincreasing. diffsnorms(A,V,D) outputs an estimate of the
%     spectral norm of A-VDV', which should be close to the machine
%     precision.
%
%
%   Reference:
%   Nathan Halko, Per-Gunnar Martinsson, and Joel Tropp,
%   Finding structure with randomness: Stochastic algorithms
%   for constructing approximate matrix decompositions,
%   arXiv:0909.4061 [math.NA; math.PR], 2009
%   (available at http://arxiv.org).
%
%
%   See also DIFFSNORMS, EIGENN, PCA, EIG, EIGS.
%

%   Copyright 2014 Mark Tygert.



%
% Check the number of inputs.
%
if(nargin < 1)
  error('MATLAB:eigens:TooFewIn',...
        'There must be at least 1 input.')
end

if(nargin > 4)
  error('MATLAB:eigens:TooManyIn',...
        'There must be at most 4 inputs.')
end

%
% Check the number of outputs.
%
if(nargout > 2)
  error('MATLAB:eigens:TooManyOut',...
        'There must be at most 2 outputs.')
end

%
% Set the inputs k, its, and l to default values, if necessary.
%
if(nargin == 1)
  k = 6;
  its = 4;
  l = k+2;
end

if(nargin == 2)
  its = 4;
  l = k+2;
end

if(nargin == 3)
  l = k+2;
end

%
% Check the first input argument.
%
if(~isfloat(A))
  error('MATLAB:eigens:In1NotFloat',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
  error('MATLAB:eigens:In1Empty',...
        'Input 1 must not be empty.')
end

%
% Retrieve and check the dimensions of A.
%
[m n] = size(A);

if(m ~= n)
  error('MATLAB:eigens:In1NotSymm',...
        'Input 1 must be square.')
end

%
% Check the remaining input arguments.
%
if(size(k,1) ~= 1 || size(k,2) ~= 1)
  error('MATLAB:eigens:In2Not1x1',...
        'Input 2 must be a scalar.')
end

if(size(its,1) ~= 1 || size(its,2) ~= 1)
  error('MATLAB:eigens:In3Not1x1',...
        'Input 3 must be a scalar.')
end

if(size(l,1) ~= 1 || size(l,2) ~= 1)
  error('MATLAB:eigens:In4Not1x1',...
        'Input 4 must be a scalar.')
end

if(~(k > 0))
  error('MATLAB:eigens:In2NonPos',...
        'Input 2 must be > 0.')
end

if((nargin > 1) && (k > n))
  error('MATLAB:eigens:In2TooBig',...
        'Input 2 must be <= the dimension of Input 1.')
end

if((nargin == 1) && (k > n))
  error('MATLAB:eigens:InTooTiny',...
        'The dimension of the input must be >= 6.')
end

if(~(its >= 0))
  error('MATLAB:eigens:In3Bad',...
        'Input 3 must be >= 0.')
end

if(l < k)
  error('MATLAB:eigens:In4ltIn2',...
        'Input 4 must be >= Input 2.')
end

%
% Check whether A is self-adjoint to nearly the machine precision
% and warn the user if not (but do NOT report an error).
%
x = rand(n,1);
y = A*x;
z = (x'*A)';

if((norm(y-z) > .1d-11*norm(y)) || (norm(y-z) > .1d-11*norm(z)))
  warning('The input matrix is not exactly self-adjoint.')
end



%
% Eigendecompose A directly if l >= n/1.25.
%
if(l >= n/1.25)

  if(~issparse(A))
    [V,D] = eig(A);
  end

  if(issparse(A))
    [V,D] = eig(full(A));
  end

%
% Rearrange the decomposition so that the absolute values
% of the eigenvalues are nonincreasing.
%
  [E,P] = sort(abs(diag(D)),'descend');
  D = diag(D);
  D = diag(real(D(P)));
  V = V(:,P);

  clear E P;

%
% Retain only the leftmost k columns of V and
% the uppermost leftmost k x k block of D.
%
  V = V(:,1:k);
  D = D(1:k,1:k);

%
% Fill the output array.
%
  if((nargout == 0) || (nargout == 1))
    varargout{1} = diag(D);
  end

  if(nargout == 2)
    varargout{1} = V;
    varargout{2} = D;
  end

  return

end


%
% Apply A to a random matrix, obtaining Q.
%
if(isreal(A))
  R = 2*rand(n,l)-ones(n,l);
end

if(~isreal(A))
  R = (2*rand(n,l)-ones(n,l)) + i*(2*rand(n,l)-ones(n,l));
end

Q = A*R;

%
% Form a matrix Q whose columns constitute a well-conditioned basis
% for the columns of the earlier Q.
%
if(its == 0)
  [Q,R,E] = qr(Q,0);
end

if(its > 0)
  [Q,R] = lu(Q);
end

%
% Conduct normalized power iterations.
%
for it = 1:its

  Q = A*Q;

  if(it < its)
    [Q,R] = lu(Q);
  end

  if(it == its)
    [Q,R,E] = qr(Q,0);
  end

end

clear E;

%
% Eigendecompose Q'*A*Q to obtain approximations to the eigenvalues
% and eigenvectors of A.
%
R = Q'*A*Q;
R = (R+R')/2;
[V,D] = eig(R);
V = Q*V;

clear Q R;

%
% Rearrange the decomposition so that the absolute values
% of the eigenvalues are nonincreasing.
%
[E,P] = sort(abs(diag(D)),'descend');
D = diag(D);
D = diag(real(D(P)));
V = V(:,P);

clear E P;

%
% Retain only the leftmost k columns of V and
% the uppermost leftmost k x k block of D.
%
V = V(:,1:k);
D = D(1:k,1:k);

%
% Fill the output array.
%
if((nargout == 0) || (nargout == 1))
  varargout{1} = diag(D);
end

if(nargout == 2)
  varargout{1} = V;
  varargout{2} = D;
end
