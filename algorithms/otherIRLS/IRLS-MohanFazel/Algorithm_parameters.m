%% -------------- PARAMETERS OF THE (s)IRLS ALGORITHM ---------------- %%



function [sr,p,rmax,fr,eta,niter,svditer,incr,gam0,gammin,tol] = Algorithm_parameters(n,r,non_zero,type)

if type == 2
    
    sr = non_zero/(n*n);
    p = non_zero; % # Measurements
    rmax = ceil(n*(1 - sqrt(1 - sr)));
    fr = r*(2*n - r)/p;
    if(fr < 0.4)
        eta = 1.1;
        niter = 500;
        svditer = 10;
        incr = 50;
    else
        eta = 1.03;
        niter = 5000;
        svditer = 50;
        incr = 100;
    end;
end;
    
gam0 = 1e-2; gammin = 1e-10; %Choose the gamma parameters - Initial and final.
tol = 1e-3; % Tolerance for convergence