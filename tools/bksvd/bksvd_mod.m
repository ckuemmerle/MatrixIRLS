function [U,S,V] = bksvd_mod(A, k, iter, bsize, center)
%--------------------------------------------------------------------------
% Randomized block Krylov iteration for truncated singular value decomposition
% Computes approximate top singular vectors and corresponding values
% Described in Musco, Musco, 2015 (http://arxiv.org/abs/1504.05477).
% =========================================================================
% Original code available at: https://github.com/cpmusco/bksvd 
% Modified by Christian Kuemmerle to include the application of function
% handles. October 15, 2017
% =========================================================================
%
% usage : 
%
%  input:
%  * A : matrix to decompose (or (1 x 4) cell where:
%                   A{1} is a function handle for the multiplication of the
%                   matrix with a vector A*v,
%                   A{2} is a function handle for the multiplication of the
%                   conjugate transpose.
%                   A{3} is the first dimension of the matrix A,
%                   A{4} is the second dimension of the matrix A.
%  * k : number of singular vectors to compute, default = 6
%  * iter : number of iterations, default = 3
%  * bsize : block size, must be >= k, default = k
%  * center : set to true if A's rows should be mean centered before the
%  singular value decomposition (e.g. when performing principal component 
%  analysis), default = false
%
%
%  output:
%  k singular vector/value pairs. 
%  * U : a matrix whose columns are approximate top left singular vectors for A
%  * S : a diagonal matrix whose entries are A's approximate top singular values
%  * V : a matrix whose columns are approximate top right singular vectors for A
%
%  U*S*V' is a near optimal rank-k approximation for A
%--------------------------------------------------------------------------

% Check input arguments and set defaults.
if nargin > 5
    error('bksvd:TooManyInputs','requires at most 5 input arguments');
end
if nargin < 1
    error('bksvd:TooFewInputs','requires at least 1 input argument');
end
if nargin < 2
    k = 6;
end
if ~iscell(A)
    k = min(k,min(size(A)));
else
    k = min(k,min(A{3},A{4}));
end

if nargin < 3
    iter = 3;
end
if nargin < 4
    bsize = k;
end
if nargin < 5
    center = false;
end
if(k < 1 || iter < 1 || bsize < k)
    error('bksvd:BadInput','one or more inputs outside required range');
end

% Calculate row mean if rows should be centered.
if ~iscell(A)
    u = zeros(1,size(A,2));
    if(center)
        u = mean(A);
    end
    l = ones(size(A,1),1);
else
   AA=A{1};
   At=A{2};
   n1=A{3};
   n2=A{4};
   u = zeros(1,n2);
   if(center)
       err('center does not work if matrix is given as function handle.');
   end
   l = ones(n1,1);
end


% We want to iterate on the smaller dimension of A.
if ~iscell(A)
    [n, ind] = min(size(A));
else
    [n, ind] = min([n1,n2]);
end
tpose = false;
if(ind == 1) 
    tpose = true;
    l = u'; 
    if ~iscell(A)
        u = ones(1,size(A,1));
        A = A';
    else
        u = ones(1,n1);
        At=A{1};
        AA=A{2};
        n2=A{3};
        n1=A{4};
    end
end

% Allocate space for Krylov subspace.
if ~iscell(A)
    K = zeros(size(A,2),bsize*iter);
    % Random block initialization.
    block = randn(size(A,2),bsize);
else
    K = zeros(n2,bsize*iter);
    % Random block initialization.
    block = randn(n2,bsize);
end
[block,R] = qr(block,0);

% Preallocate space for temporary products.
if ~iscell(A)
    T = zeros(size(A,2),bsize);
else
    T = zeros(n2,bsize);
end

% Construct and orthonormalize Krlov Subspace. 
% Orthogonalize at each step using economy size QR decomposition.
if ~iscell(A)
    for i=1:iter
        T = A*block - l*(u*block);
        block = A'*T - u'*(l'*T);
        [block,R] = qr(block,0);
        K(:,(i-1)*bsize+1:i*bsize) = block;
    end
else    
    for i=1:iter
        T = feval(AA,block) - l*(u*block);
        block = feval(At,T) - u'*(l'*T);
        [block,R] = qr(block,0);
        K(:,(i-1)*bsize+1:i*bsize) = block;
    end
end
[Q,R] = qr(K,0);

% Rayleigh-Ritz postprocessing with economy size dense SVD.
if ~iscell(A)
    T = A*Q - l*(u*Q);
else
    T = feval(AA,Q) - l*(u*Q);
end

[Ut,St,Vt] = svd(T,0);
S = St(1:k,1:k);
if(~tpose)
    U = Ut(:,1:k);
    V = Q*Vt(:,1:k);
else
    V = Ut(:,1:k);
    U = Q*Vt(:,1:k);
end

end