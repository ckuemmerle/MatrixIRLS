function times = computeS(n_1,n_2,r,m,s)

t1 = zeros(s,1); t2 = zeros(s,1);
for i = 1:s
    X = orth(randn(n_1,r));
    Y = orth(randn(n_2,r));
    M = prob_matrix(n_1,n_2,r,m);
    M_E = M.sparse;
    E = spones(M_E);
    [t1(i),t2(i)] = compS(X,Y,M_E,E);
%     fprintf('LS solution without projection: %3.2f s \n',t1);
%     fprintf('LS solution with projection: %3.2f s \n',t2);
end
times(1) = min(t1);
times(2) = min(t2);



function [time1,time2] = compS(X,Y,M_E,E)

tau = tic;
r = size(X,2);
C = X' * ( M_E ) * Y; 
C = C(:); %right-hand side vector
time0 = toc(tau);

%original part
t = tic;
for i = 1:r
    
    for j = 1:r
                
                ind = (j-1)*r + i ;
                temp = X' * (  (X(:,i) * Y(:,j)').*E ) * Y ;
                A(:,ind) = temp(:) ;
				                
    end
    
end

S = A\C ; 
out1 = reshape(S,r,r) ;
time1 = toc(t) + time0;

%improved part
t = tic;
B = E;
[I,J] = find(E);
l = nnz(E);
A = zeros(r^2,r^2);
for i = 1:r
    
    for j = 1:r
                
                ind = (j-1)*r + i ;
                val = partXY(X(:,i)',Y(:,j)',I,J,l);
                updateSval(B,val,l);
                %val = sparse(I,J,val',n,m);
                temp = X' * B * Y;
                A(:,ind) = temp(:) ;
				
    end
    
end


S = A\C ; 
out2 = reshape(S,r,r) ;
time2 = toc(t) + time0;

% %improved part II
% t = tic;
% [n_1,~] = size(X);
% [n_2,r] = size(Y);
% omega = find(E);
% m = nnz(E);
% b = M_E(omega);
% 
% for i = 1:m
%   
%   k = ceil(omega(i)/(n_1));
%   l = mod(omega(i)-1,n_2)+1;
%     
%   temp = X(l,:)*Y(k,:)';
%   A(i,:) = temp(:);
%     
% end
% 
% 
% S = A\b; 
% out3 = reshape(S,r,r);
% time3 = toc(t)
% 
% norm(out3-out1,'fro')/norm(out3,'fro')
