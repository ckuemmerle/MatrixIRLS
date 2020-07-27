%%------------ GRADIENT PROJECTION ---------------------%%
% ----- KARTHIK MOHAN (karna@uw.edu)---------%
%------------ LAST UPDATE 8/28/2012 -------------------%

function[Xnew,err,terr,k] =   grad_proj(B,L,Xprev,V,D1,m,n,alpt,betat,kmax)

%% PARAMETERS

%Initialize
step = 1;

Xnew = Xprev;


tol = 1e-4; %TOLERANCE

k = 1; err = 10; 
terr = zeros(kmax,1);


while(k < kmax && err > tol)

Y = Xnew;
Xold = Xnew;

temp = Y - (2*step)/L*((Y*V)*(D1*V') + Y); %Renormalized W being used here
%Note - Step size here is actually step/L - The 2 comes out of the gradient.

Xnew = temp - mcop(temp,m,n,alpt,betat) + B; %GRADIENT PROJECTION
err = norm(Xnew - Xold,'fro')/norm(Xold,'fro');

terr(k,1) = err;
k = k + 1;
end;

return;
