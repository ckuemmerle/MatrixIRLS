function [x,histout,fail,Xout,itc,time,msg,debug] = LRGeomCG(prob, opts, x0)
% LRGEOMCG    Low-rank matrix completion by Geometric CG.
%   Uses Armijo rule, polynomial linesearch on embedded submanifold of
%   fixed rank matrices.
%
% Input: prob     = problem instance, see MAKE_PROB.
%        opts     = options, see DEFAULT_OPTS.
%        x0       = starting guess.
%
% Output: x       = solution.
%         histout = iteration history. Each row of histout is
%                   [rel_norm(grad), rel_err_on_omega, relative_change, ...
%                        number of step length reductions, restarts]
%         fail    = flag for failure
%
% See also SMALL_EXAMPLE, MAKE_PROB, DEFAULT_OPTS.

% (C) Bart Vandereycken, 2011-2012.
% Adapted from steep by C. T. Kelley, Dec 20, 1996.
%
% Paper: Low-rank matrix completion by Riemannian optimization.
%
% Info and bugs: bart.vandereycken@epfl.ch
% =========================================================================
% Code modified by Christian Kuemmerle:
% Version of LRGeomCG to serve as inner method in LRGeomCG_pursuit.
% Adapted to interface with desired format for algorithmic comparisons.
% Feburary 2021.
% =========================================================================


%beta_type = 'F-R';
beta_type = 'P-R';

fail = true;

norm_M_Omega = norm(prob.data);

ORTH_VALUE = 0.1; % the search directions should be almost orthogonal


itc = 1; xc = x0;
fc = F(prob,xc);
gc = grad(prob,xc);
ip_gc = ip(xc,gc,gc);
rel_grad = sqrt(ip_gc);
fold = 2*fc; reschg = -1;
% first search-dir is steepest gradient
dir = scaleTxM(gc,-1);

opts.rel_grad_tol = max(opts.rel_grad_tol, rel_grad*opts.rel_grad_decrease_factor);


ithist=zeros(opts.maxit,4);
beta = 0;

c1 = 0.001;
c2 = 0.1;
prtlevel = opts.verbosity;

msg = 'Max iterations reached.'; 
ithist(itc,1) = rel_grad;
ithist(itc,2) = sqrt(2*fc)/norm_M_Omega;
ithist(itc,3) = reschg;
ithist(itc,4) = 0;
if prob.has_testing
    on_gamma = partXY(xc.S'*xc.U',xc.V', prob.Gamma_i, prob.Gamma_j, prob.m);
    ithist(itc,5) = 0.5*norm( on_gamma - prob.test_data' ) / norm(prob.test_data);
else
    ithist(itc,5) = Inf;
end
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    Xout=cell(1,opts.maxit);% Xout=zeros(prob.n1,prob.n2,opts.maxit);
end
time=zeros(1,opts.maxit);
debug{1} = xc;

tic
for itc=2:opts.maxit
    
    
    % perform line search with initial guess from linearized problem
    tinit = exact_search_onlyTxM(prob, xc, dir);
    
    fc_old = fc;
    if opts.strong_wolfe
        [lambda, xc_new, fc, gc_new, succ, iarm] = linesch_sw(prob, xc, fc, gc, dir, c1, c2, tinit, prtlevel);
    else
        [xc_new,fc,gc_new, succ,~,iarm,lambda] = armijo_search(prob, xc,fc,gc,dir, tinit);
    end
    
    
    % line search failed, reset to  steepest descent and try again
    if ~succ && beta ~= 0
        if opts.verbosity > 0; disp('Line search failed on CG. Resetting to gradient.'); end
        beta = 0;
        % if we did cg, reset to steepest descent and try again
        dir = scaleTxM(gc,-1);
        tinit = exact_search_onlyTxM(prob, xc, dir);
        if opts.strong_wolfe
            [~, xc_new, fc, gc_new, succ, iarm] = linesch_sw(prob, xc, fc_old, gc, dir, c1, c2, tinit, prtlevel);
        else
            [xc_new,fc,gc_new, succ,~,iarm] = armijo_search(prob, xc,fc_old,gc,dir, tinit);
        end
    end
    
    debug{itc} = xc_new;
    
    % if it still fails (beta is always 0 here) then we give up
    if ~succ
        x = xc_new;
        ithist(itc,1) = rel_grad;
        ithist(itc,2) = sqrt(2*fc)/norm_M_Omega;
        ithist(itc,3) = reschg;
        ithist(itc,4) = iarm;
        if prob.has_testing
            on_gamma = partXY(x.S'*x.U',x.V', prob.Gamma_i, prob.Gamma_j, prob.m);
            ithist(itc,5) = 0.5*norm( on_gamma - prob.test_data' ) / norm(prob.test_data);
        else
            ithist(itc,5) = Inf;
        end
        histout=ithist(1:itc,:);
        msg = 'Line search failed.';
        if opts.verbosity > 0; disp(msg); end
        if isfield(opts,'saveiterates') && opts.saveiterates == 1
            x.Xout=Xout(1:itc-2);
            Xout= x;
        else
            x.Xout = {{xc_new.U*xc_new.S,xc_new.V}};
            Xout= x;
        end
        itc=itc-2;
        time=time(1:itc);
        return
    end
    if isfield(opts,'saveiterates') && opts.saveiterates == 1
        Xout{itc-1}= {xc_new.U*xc_new.S,xc_new.V};
    end
    time(itc-1) = toc;
    % grad(new x)
    %gc_new = grad(prob,xc_new);
    ip_gc_new = ip(xc_new,gc_new,gc_new);
    
    % Test for convergence
    if sqrt(2*fc) < opts.abs_f_tol
        msg = 'Abs_f_tol reached.';
        if opts.verbosity > 0; disp(msg); end
        fail = false;
        break;
    end
    if sqrt(2*fc)/norm_M_Omega < opts.rel_f_tol
        msg = 'Rel_f_tol reached.';
        if opts.verbosity > 0; disp(msg); end
        fail = false;
        break;
    end
    
    rel_grad = sqrt(ip_gc_new);
    if rel_grad < opts.rel_grad_tol
        msg = 'Rel_grad_tol reached.';
        if opts.verbosity > 0; disp(msg); end
        fail = false;
        break;
    end
    
    % for 'stagnation stopping criterion' 
    reschg = abs(1-sqrt(2*fc)/sqrt(2*fold) );  % LMARank's detection
    %reschg = abs(sqrt(2*fc) - sqrt(2*fold)) / max(1,norm_M_Omega);
    if itc > 4 && reschg < opts.rel_tol_change_res
        msg = 'Rel_tol_change_res reached.';
        if opts.verbosity > 0; disp(msg); end
        fail = true;
        break;
    end
    
    
    if itc > 10
        R1 = triu(qr([xc_new.U*xc_new.S xc.U*xc.S],0));
        R2 = triu(qr([xc_new.V -xc.V],0));
        
        %rel_change_x = norm(xc_new.U*diag(xc_new.sigma)*xc_new.V' - xc.U*diag(xc.sigma)*xc.V', 'fro') / ...
        %  norm(xc_new.U*diag(xc_new.sigma)*xc_new.V', 'fro');
        rel_change_x = norm(R1*R2','fro')/norm(diag(xc_new.S),2);
        
        if rel_change_x < opts.rel_tol_change_x
            msg = 'Rel_tol_change_x reached.';
            if opts.verbosity > 0; disp(msg); end
            fail = true;
            break;
        end
    end
    
    % new dir = -grad(new x)
    %           + beta * vectTransp(old x, old dir, tmin*old dir)
    
    % differentiated retraction, makes no difference (it equal up to
    % second-order to projected)
    if opts.differentiated_transport
        [xc_new2, gc_old_2] = transpVect_differentiated(prob,xc,gc, scaleTxM(dir,lambda),0);  % the transported will be in another basis... so project again
        gc_old = transpVect(prob,xc_new2,gc_old_2,xc_new,1);
        [xc_new2, dir_2] = transpVect_differentiated(prob,xc,dir, scaleTxM(dir,lambda),0);
        dir = transpVect(prob,xc_new2,dir_2,xc_new,1);
    else
        gc_old = transpVect(prob,xc,gc,xc_new,1);
        dir = transpVect(prob,xc,dir,xc_new,1);
    end
%     
    
    
    % we test how orthogonal the previous gradient is with
    % the current gradient, and possibly reset the to the gradient
    ip_gc_old_new = ip(xc_new,gc_old,gc_new);
    orth_grads = ip_gc_old_new / ip_gc_new;
    
    if abs(orth_grads) >= ORTH_VALUE 
        if opts.verbosity > 1; disp('New gradient is almost orthogonal to current gradient. This is good, so we can reset to steepest descent.'); end
        beta = 0;
        dir = plusTxM(gc_new, dir, -1, beta);
        
    else % Compute the CG modification
        % beta
        if strcmp(beta_type, 'F-R')  % Fletcher-Reeves
            beta = ip_gc_new / ip_gc;
            % new dir
            dir = plusTxM(gc_new, dir, -1, beta);
        elseif strcmp(beta_type, 'P-R')  % Polak-Ribiere+
            % vector grad(new) - transported grad(current)
            beta = (ip_gc_new - ip_gc_old_new) / ip_gc;
            beta = max(0,beta);
            dir = plusTxM(gc_new, dir, -1, beta);
        end
    end
    
    % check if dir is descent, if not take -gradient (i.e. beta=0)
    g0 = ip(xc_new,gc_new,dir);
    if g0>=0
        if opts.verbosity > 1
            disp('New search direction not a descent direction. Resetting to -grad.');
        end
        dir = scaleTxM(gc_new, -1);
        beta = 0;
    end
    
    
    % update _new to current
    gc = gc_new;
    ip_gc = ip_gc_new;
    xc = xc_new;
    fold = fc;
    
    ithist(itc,1) = rel_grad;
    ithist(itc,2) = sqrt(2*fc)/norm_M_Omega;
    ithist(itc,3) = reschg;
    ithist(itc,4) = iarm;
    if prob.has_testing
        on_gamma = partXY(xc.S'*xc.U',xc.V', prob.Gamma_i, prob.Gamma_j, prob.m);
        ithist(itc,5) = 0.5*norm( on_gamma - prob.test_data' ) / norm(prob.test_data);
    else
        ithist(itc,5) = Inf;
    end
end
itc=itc-1;
time=time(1:itc);
x = xc_new;
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    x.Xout=Xout(1:itc);
    Xout= x;
else
    x.Xout = {{xc_new.U*xc_new.S,xc_new.V}};
    Xout= x;
end

ithist(itc,1) = rel_grad;
ithist(itc,2) = sqrt(2*fc)/norm_M_Omega;
ithist(itc,3) = reschg;
ithist(itc,4) = iarm;
if prob.has_testing
    on_gamma = partXY(x.S'*x.U',x.V', prob.Gamma_i, prob.Gamma_j, prob.m);
    ithist(itc,5) = 0.5*norm( on_gamma - prob.test_data' ) / norm(prob.test_data);
else
    ithist(itc,5) = Inf;
end
histout=ithist(1:itc,:);
       
if opts.verbosity > 0; disp(msg); end
end

