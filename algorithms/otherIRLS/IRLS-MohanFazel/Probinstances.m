%% ------- PROBLEM INSTANCES TO TEST IRLS-q, sIRLS-q ----------- %%

% ----- LAST UPDATE: 8/28/2012 ----------------------------- %

function [m,n,sr,p,r,rmax,fr,eta,svditer,incr,niter] = Probinstances(h,j);

etamin = 1.02; etamax = 1.15;
frmin = 0.17; frmax = 0.87;

if(h == 1)
     
%% Easy problem instances FR < 0.4
niter = 500; %Max # outer iterations   
svditer = 10; 
incr = 50;  %Randomized SVD parameter to switch between QR and svds 

if(j == 1)
    m = 100; n = 100; sr = 0.57; p = round(m*n*sr); r = 10; 
    eta = 1.1;
end;
if(j == 2)
    m = 200; n = 200; sr = 0.39; p = round(m*n*sr); r = 10; 
    eta = 1.1;
end;
if(j == 3)
    m = 500; n = 500; sr = 0.2; p = round(m*n*sr); r = 10; 
    eta = 1.1;
end;
if(j == 4)
    m = 500; n = 500; sr = 0.12; p = round(m*n*sr); r = 10; 
    eta = 1.1;
end;
if(j == 5)
    m = 1000; n = 1000; sr = 0.12; p = round(m*n*sr); r = 10; 
    eta = 1.1;
end;
if(j == 6)
    m = 1000; n = 1000; sr = 0.39; p = round(m*n*sr); r = 50; 
    eta = 1.1;
end;
if(j == 7)
    m = 1000; n = 1000; sr = 0.12; p = round(m*n*sr); r = 20; 
    eta = 1.1;
end;
if(j == 8)
    m = 2000; n = 2000; sr = 0.12; p = round(m*n*sr); r = 20; 
    eta = 1.1;
end;
if(j == 9)
    m = 2000; n = 2000; sr = 0.12; p = round(m*n*sr); r = 40; 
    eta = 1.1;    
end;
rmax = ceil(n*(1 - sqrt(1 - sr)));
else
 
%% Hard problem instances FR > 0.4
      
niter = 10000; %Max # outer iterations 
svditer = 50;  %# initial iterations for which  greater than rank r approx
              %of X^k is made.
incr = 100;   %Randomized SVD parameter to switch between QR and svds 

if(j == 1)
    m = 40; n = 40; sr = 0.5; p = round(m*n*sr); r = 9; 
    svditer = 400; eta = 1.03; incr = 4000;
end;

if(j == 2)
    m = 100; n = 100; sr = 0.3; p = round(m*n*sr); r = 14; 
    svditer = 400; eta = 1.03; incr = 4000;
end;
if(j == 3)
    m = 500; n = 500; sr = 0.1; p = round(m*n*sr); r = 20; 
    eta = 1.03;
end;

if(j == 4)
    m = 1000; n = 1000; sr = 0.1; p = round(m*n*sr); r = 20; 
    eta = 1.03;
end;
if(j == 5)
    m = 1000; n = 1000; sr = 0.08; p = round(m*n*sr); r = 20; 
    eta = 1.03;
end;
if(j == 6)
    m = 1000; n = 1000; sr = 0.07; p = round(m*n*sr); r = 20; 
    eta = 1.03;
   
end;
if(j == 7)
    m = 1000; n = 1000; sr = 0.06; p = round(m*n*sr); r = 20; 
    eta = 1.03;
end;

if(j == 8)
    m = 1000; n = 1000; sr = 0.1; p = round(m*n*sr); r = 30; 
    eta = 1.03;
end;
if(j == 9)
    m = 1000; n= 1000; sr = 0.2; p = round(m*n*sr); r = 50; 
    eta = 1.03;
end;


rmax = ceil(n*(1 - sqrt(1 - sr)));
end;
fr = r*(2*n - r)/p;
%eta = etamin + (etamax - etamin)/(frmin - frmax)*(fr - frmax);
%fr = r*(2*n - r)/p;
return