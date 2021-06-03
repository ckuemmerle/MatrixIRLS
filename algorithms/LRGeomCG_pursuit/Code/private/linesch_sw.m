function [alpha, x, f, grad, succ, nsteps] = ...
           linesch_sw(prob, x0, f0, grad0, d, c1, c2, tinit, prtlevel)
%  based on Nocedal-Wright line search, see their book
%  strong Wolfe line search with cubic interpolation
%  not recommended for BFGS, especially if function may be nonsmooth: 
%    use the far simpler and preferable linesch_ww instead
%  recommended for use with CG, where strong Wolfe condition needed for
%  convergence analysis
% call:  [alpha, x, f, grad, fail, nsteps] = ...
%          linesch_sw(x0, f0, grad0, d, pars, c1, c2, prtlevel);
%  written by M. Overton (overton@cs.nyu.edu)
%  input
%   x0:      intial point
%   f0:      function value at x0
%   g0:      gradient at x0
%   d:       search direction (called p in Nocedal and Wright's book) 
%   c1: Wolfe parameter for the condition 
%          f(x0 + alpha d) < = f0 + c1*alpha*grad0'*d   (default 0)
%        in theory, should be positive, but 0 works fine in practice
%   c2: Wolfe parameter for the (strong) condition 
%          abs((grad f)(x0 + alpha d)'*d) < = -c2*grad0'*d
%        these should satisfy 0 <= c1 <= c2 <= 1     (default 0.5)
%        usually, c1 should be small but c2 should not, but set c2 small
%         if an accurate line search is desired
%        for BFGS, only need c2 < 1 for update to be positive definite 
%        for nonlinear conjugate gradient, need c2 < 0.5 for the 
%         Fletcher-Reeves convergence analysis to apply
%   prtlevel: 2 for info when line search fails (default), 0 for no printing
%      
%  output:
%   alpha: steplength satisfying Wolfe conditions
%   x:     x0 + alpha*d
%   f:     f(x0 + alpha d)
%   grad:  (grad f)(x0 + alpha d)
%   fail:  0 for normal termination (Wolfe conditions satisfied)
%          1 if too many steps were taken in lszoom (Wolfe conds not satisfied)
%          -1 if minimizer was not bracketed (too many expansion steps, 
%                                              function may be unbounded below)
%   nsteps: number of steps taken in lszoom
% Adapted from linesch_sw by Michael Overton, GPL license.


if nargin < 6  % check if the optional Wolfe parameters were passed
    c1 = 0.01;  % should be positive in theory, not necessary in practice
end
if nargin < 7
    c2 = 0.49;  % should be < 1/2 in theory
end
if nargin < 9
    prtlevel = 2;
end
if c1 < 0 || c1 > c2 || c2 > 1
   warning('linesch_sw: Wolfe parameters are supposed to satisfy 0 <= c1 <= c2 <= 1')
end
g0 = ip(x0,grad0,d);  % g0 = grad0'*d;
if g0 >= 0
   error('linesch_sw: g0 must be negative, not a descent direction')
end
old = 0; fold = f0; gold = g0;
new = tinit;
nexpand = max([50  -round(log2(sqrt(ip(x0, d, d))))]);

store = [];

for k = 1:nexpand
    [fnew, xnew, store, R_gradnew, gnew] = F_line_search(prob,x0,d,new,store);
   %xnew = x0 + new*d' %[fnew,gradnew] = feval(fgname, xnew, pars); gnew = gradnew'*d;
   if fnew > f0 + c1*new*g0 || ((fnew >= fold) && k > 1)   % gone too far
      [alpha, x, f, grad, fail, nsteps] = lszoom(prob, store, old, new, ...
          fold, fnew, gold, gnew, f0, g0, x0, d, c1, c2, prtlevel);
      succ = ~fail;
      return
   end
   if abs(gnew) <= -c2*g0   % quit, conditions are satisfied first time
      alpha = new; x = xnew; f = fnew; grad = R_gradnew; fail = false; nsteps = k;
      succ = ~fail;
      return
   end
   if gnew >= 0     % gone past minimizer, but new point is better than old
      [alpha, x, f, grad, fail, nsteps] = lszoom(prob, store, new, old, ... %note reversal
          fnew, fold, gnew, gold, f0, g0, x0, d, c1, c2, prtlevel);
      succ = ~fail;
      return
   end
   old = new;       % minimizer not bracketed yet
   fold = fnew;
   gold = gnew;
   new = 2*new;     
end % loop
if prtlevel > 1
    fprintf('linesch_sw: minimizer was not bracketed in %d expansion steps, ', nexpand)
    fprintf('function may be unbounded below\n')
end
alpha = new; 
x = xnew; 
f = fnew; 
grad = gradnew; 
succ = false;
nsteps = 0;  % lszoom was never called