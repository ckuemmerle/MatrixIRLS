%% ---- Test file for IRLS-q ---------------- %%

%% ----- This is the code associated with  the paper:
% ----- "Iterative Reweighted Algorithms for Matrix Rank Minimization"
% ----- Karthik Mohan (karna@uw.edu) and Maryam Fazel (mfazel@uw.edu).

% ----- LAST UPDATE: 8/28/2012 --------------%

clc; clear; close all;

%% PROBLEM SETUP

Problem_setup; % Set the problem up.
nrg = 1; % Enter the number of random generations of the data required to average the results.


%% CHECK IF PROBLEM PARAMETERS ARE IN CORRRECT RANGE
while(type < 1 || type > 2 || abs(type - floor(type)) > 0)
    fprintf('Choose 1 for testing synthetic data and 2 for testing real data');
    while(type < 1 || type > 2 || abs(type - floor(type)) > 0)
    type = input('\n Enter either 1 or 2 ');
    end;
    ins = 0;
    while(ins < 1 || ins > 9 || abs(ins - floor(ins)) > 0)
    ins = input('Choose a problem instance between 1 and 9: ');
    end;
    h = 0;
    while h < 1 || h > 2 || abs(h - floor(h)) > 0
    fprintf('\n Choose between easy and hard problem instances.\n');
    h = input('\n Enter an integer between 1 and 2: ');
    end;
end;

if type == 1
while(ins < 1 || ins > 9 || abs(ins - floor(ins)) > 0)
   fprintf('\n The problem instance has to be an integer between 1 and 9.\n');
   ins = input('\n Enter an integer between 1 and 9:  ');
end;

while h < 1 || h > 2 || abs(h - floor(h)) > 0
    fprintf('\n Choose between easy and hard problem instances.\n');
    h = input('\n Enter an integer between 1 and 2: ');
end;

end;

if type == 2
if(exist('M.mat') == 0)
    return;
end;
end;

while(nrg <= 0 || abs(nrg - floor(nrg)) > 0)
    fprintf('\n The number of random generations has to be positive integer\n');
    nrg = input('\n Enter an integer greater than 0:  ');
end;



%% CHOOSE ALL REMAINING PARAMETERS
if type == 1
[m,n,sr,p,r,rmax,fr,eta,svditer,incr,niter] = Probinstances(h,ins);
M = zeros(m,n);
gam0 = 1e-2; gammin = 1e-10; %Choose the gamma parameters - Initial and final.
tol = 1e-3; % Tolerance for convergence
end;

if(type == 2)
    non_zero = size(M,1);
    [sr,p,rmax,fr,eta,niter,svditer,incr,gam0,gammin,tol] = Algorithm_parameters(n,r,non_zero,type);
end;

q = 0; %sIRLS-q: Choose a q between 0 and 1

while(q < 0 || q > 1)
    q = input('\n Enter a real number between 0 and 1:  ');
end;

rknown = 1; % CHOOSE 1 if the Algorithm is allowed to use the
            % information on the rank of the true solution.
            % CHOOSE 2 if the Algorithm is unware of the 
            % rank of the true solution.
while(rknown <0 || rknown > 1 || abs(rknown - floor(rknown)) > 0)
    rknown = input('\n Enter either 0 or 1:  ');
end;


fprintf('\n -------------------');
fprintf('\n Algorithm begins...');
fprintf('\n -------------------\n\n');
[NS, avgerr,avgiterno, TT,timeperiter, TTcpu,Xalgo] = irls_q(m,n,sr,r,rmax,rknown,eta,gam0,gammin,q,tol,nrg,niter,svditer,incr,type,M);

%% ----------- OUTPUT ---------------------- %%

if(type == 1)
fprintf('\n\n m = %d, n = %d, r = %d, p = %d, samp.ratio = %3.2f, freedom = %3.2f, eta = %1.3f \n', m,n,r,p,sr,fr,eta);
fprintf(' NS = %d, Avg Rec Err. = %0.5f,Avg # Iters = %d, Avg clock time = %3.2f,Clock time/iter = %3.3f Avg cpu time = %3.2f \n\n\n', NS, avgerr,avgiterno, TT,timeperiter, TTcpu);
else
fprintf('\n\n m = %d, n = %d, r = %d, p = %d, samp.ratio = %3.2f, freedom = %3.2f, eta = %1.3f \n', m,n,r,p,sr,fr,eta);
fprintf(' # Iters = %d, Clock time = %3.2f, Clock time/iter = %3.3f Cpu time = %3.2f \n\n\n', avgiterno, TT,timeperiter, TTcpu);
end;


fprintf('\n The completed matrix is given by Xalgo.mat ...\n');
