%%----------- RANDOMIZED ALGORITHM FOR SVD COMPUTATION ----------------%%
%   IMPLEMENTATION BASED ON :              %
%   PAPER BY TROPP et al:                  %
%   http://arxiv.org/PS_cache/arxiv/pdf/0909/0909.4061v1.pdf %

%------------- SECOND LAST UPDATE: 9/14/2011 --------------- %
%------------- LAST UPDATE: 6/26/2012 ---------------------- %

function [U,S,V] = rand_svd(A,k,iter_no,svditer,incr);

n = size(A,2);
inc = k;


%STEP 1 - CAPTURE THE RANGE SPACE OF A
O = randn(n,k+inc);
Y0 = A*O;
Y1 = A*(A'*Y0);
%Y2 = A*(A'*Y1);
P = [Y0 Y1];

%STEP 2 - COMPUTE LEFT SINGULAR VECTORS FOR P

if(iter_no > svditer+incr)
[Q1,R] = qr(P,0);
Q = Q1(:,1:k);
else
[Q,S1,V1] = svds(P,k);
end;

%STEP 3 - COMPUTE St
St = A'*Q;

%STEP 4 - COMPUTE SVD OF S'
[U1,S,V] = svds(St',k);
% diff=k-size(S,1)
% sizeS=size(S,1)
% if diff>0
%     s = diag(S);
%     size(S,1)+1:k
%     s = diag(S);
%     s(size(S,1)+1:k)=zeros(diff,1);
%     S=diag(s);
% end
% sizeS=size(S,1)
s = diag(S);

%STEP 5 - GET LEFT SINGULAR VECTORS OF A

U = Q*U1;



