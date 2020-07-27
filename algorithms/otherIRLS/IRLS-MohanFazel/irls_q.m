%% -------------- IRLS-q (0 <= p <= 1) algorithm -------------------- %

%% ----- This is the code associated with  the paper:
% ----- "Iterative Reweighted Algorithms for Matrix Rank Minimization"
% ----- Karthik Mohan (karna@uw.edu) and Maryam Fazel (mfazel@uw.edu).

% -------------- LAST UPDATE: 8/28/2012 ------------------------------ %


function [NS, avgerr,avgiterno, TT,timeperiter, TTcpu,Xnew] = irls_q(m,n,sr,r,rmax,rknown,eta,gam0,gammin,q,tol,nrg,niter,svditer,incr,type,M)



%% PARAMETERS
 
if(rknown == 1)
countstart = r; %-- r if rank known
else
countstart = rmax; %-- rmax if rank NOT known
end;
TT = 0;TTcpu = 0; NS = 0;%# succesful instances
err = 0; avgerr = 0; avgiterno = 0; 


%% TYPE = 1 : SYNTHETIC DATA
if(type == 1)    
    
for(ng = 1: nrg)


%% GENERATE MATRIX COMPLETION OPERATOR AND LOW RANK MATRIX X0    
Y1 = randn(m,r); Y2 = randn(n,r);
X01 = Y1*Y2'; X0 = X01/norm(X01); gam = gam0*norm(X0);

OM = binornd(1,sr,m,n);
[alp,beta] = find(OM==1); %vectors defining support Omega
[betat,alpt] = find(OM' == 1); 
p1 = size(alp,1); %NUMBER OF MEASUREMNTS TAKEN ~= p*n*n
CO = size(n,1);


%% MEASUREMENTS

B = zeros(m,n);
for(i = 1:p1)
    B(alpt(i),betat(i)) = X0(alpt(i),betat(i));
end;


%% IRLS FOR MATRIX COMPLETION PROBLEM


%FIRST ITERATION OF IRLS
k = 1;
Xnew = B; 
svditer =  niter; %Parameter in rand_svd
count = countstart;
L = 2; %Initial Lipschitz constant
V = 0; D1 = 0;
extra_rank = 0;
tstart = cputime;
startclock = clock; 
fprintf('\n');
while(k<niter)
    
    [Xnew,err,terr,l] = grad_proj(B,L,Xnew,V,D1,m,n,alpt,betat,500);  
    [U,S,V] = rand_svd(Xnew,count,k,svditer,incr);
     s = diag(S);
     g = gam*ones(count,1); s = s(1:count,1);
     D1 = diag( (g.^(1 - q/2))./((s.*s + g).^(1 - q/2))- ones(count,1));
     V = V(:,1:count);
    
     count = min(size(find(s > max(s)*1e-2),1)+extra_rank,rmax); 
     % Estimating rank to truncate SVD in next iteration
     if(rknown == 1)
         count = r;
     end;
      
    L = 2; %Lipschitz constant
    err = norm(Xnew - X0,'fro')/norm(X0,'fro');
    gam = max(gam/eta,gammin);
    k = k + 1;
    if(mod(k,40) == 0)
        fprintf('.');
        if(err < tol)
            break;
        end;
    end;
end;

if(err < tol)
avgerr = err + avgerr;
avgiterno = k + avgiterno;
TTcpu = TTcpu + cputime - tstart; 
TT = TT + etime(clock,startclock);
NS = NS + 1;
end;

end;


%% TYPE = 2 : USER-INPUT DATA

else %IF TYPE = 2

B = zeros(m,n);
for(i = 1:size(M,1))
    B(M(i,1),M(i,2)) = M(i,3);
end;
alpt = M(:,1); betat = M(:,2);
gam = gam0*norm(B);

%FIRST ITERATION OF IRLS
k = 1;
Xnew = B; Xold = B;
svditer =  niter; %Parameter in rand_svd
count = countstart;
L = 2; %Initial Lipschitz constant
V = 0; D1 = 0;
extra_rank = 0;
tstart = cputime;
startclock = clock; 
fprintf('\n');
while(k<niter)
    
    [Xnew,err,terr,l] = grad_proj(B,L,Xnew,V,D1,m,n,alpt,betat,500);  
    [U,S,V] = rand_svd(Xnew,count,k,svditer,incr);
     s = diag(S);
     g = gam*ones(count,1); s = s(1:count,1);
     D1 = diag( (g.^(1 - q/2))./((s.*s + g).^(1 - q/2))- ones(count,1));
     V = V(:,1:count);
    
     count = min(size(find(s > max(s)*1e-2),1)+extra_rank,rmax); 
     % Estimating rank to truncate SVD in next iteration
     if(rknown == 1)
         count = r;
     end;
      
    L = 2; %Lipschitz constant
    err = norm(Xnew - Xold,'fro')/norm(Xnew,'fro');
    Xold = Xnew;
    gam = max(gam/eta,gammin);
    k = k + 1;
    if(mod(k,40) == 0)
        fprintf('.');
        if(err < tol)
            break
        end;
    end;
end;

TT = cputime - tstart;
TTcpu = TT;
avgerr = err;
avgiterno = k;
timeperiter = TT/avgiterno;
NS = 1;
end;

if(type == 1)
if(NS > 0)
TTcpu = TTcpu/NS;
TT = TT/NS;
avgerr = avgerr/NS;
avgiterno = avgiterno/NS;
timeperiter = TT/avgiterno;
else
TTcpu = cputime - tstart;
TT = etime(clock,startclock);
avgerr = err;
avgiterno = k;
timeperiter = TT/avgiterno;    
end;
end;

   