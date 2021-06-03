function [x,histout,stats,Xr,outs] = LRGeomCG_pursuit(prob, opts)
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
% This code implements a version of the Riemannian pursuit algorithm
% presented in the papers:
% [1] André Uschmajew, Bart Vandereycken, "Greedy rank updates combined with Riemannian descent methods 
% for low-rank optimization", In Sampling Theory and Applications (SampTA), 
% 2015 International Conference on, pages 420–424. IEEE, 2015.
% [2] 
% Code modified by Christian Kuemmerle:
% - save intermediate iterates (if applicable)
% - some adaptations to provide output to interface with other methods
% February 2021.
% =========================================================================
t = tic();
rank_increase = opts.rank_increase;

p = length(prob.Omega) / prob.n1 / prob.n2;

opts_svds.tol = 1e-10;
opts_cg = opts;
opts_cg.verbosity = 0;

% first iteration
prob.r = rank_increase;
M_omega = sparse(prob.Omega_i, prob.Omega_j, prob.data,prob.n1,prob.n2,prob.m);
if opts.rand_pca
    [U,S,V] = pca(M_omega/p,prob.r,true); flag = false;
else
    [U,S,V,flag] = svds(M_omega/p, prob.r, 'L', opts);
    if flag, warning('svds did not converge'); end
end
x0.V = V;
x0.S = S;
x0.U = U;
x0 = prepx(prob, x0);

hist = []; chg_ranks = [];
n = size(x0.U,1); m = size(x0.V,1);

M_omega_old = M_omega;
if opts.opt_S==2    
    M_omega = sparse(prob.opt_S_fast.Omega_i, prob.opt_S_fast.Omega_j, prob.opt_S_fast.data,prob.n1,prob.n2,prob.opt_S_fast.m);        
end
ranklist = [];
its = 1;
while prob.r <= opts.max_rank
    fprintf(' Outer iteration: %3i. ', its)
    [x_cg,hist_cg,fail_cg,outs_c,N_c,tm_c,msg_cg] = LRGeomCG1(prob,opts_cg,x0);
    fprintf('Finished (%28s) Iterations = %3i, Relative error = %d, Relative gradient = %d \n', msg_cg, size(hist_cg,1), hist_cg(end,2), hist_cg(end,1))
    hist = [hist; hist_cg];
    chg_ranks = [chg_ranks size(hist,1)];
    
    if prob.r + rank_increase > opts.max_rank
        break
    end
    
%     if its == 12
%         rank_increase = rank_increase+1
%     end
%     
    %f1 = F(prob, x_cg)
    
    % increase rank and make a search direction
    if opts.increase_with_residual
        d = x_cg.on_omega - prob.data;
        updateSval(prob.temp_omega, d, length(d));
        grad_fullspace = prob.temp_omega;
         % todo orth proj: does not make a difference in practice. If
         % previous problem is solved to zero gradient, it also does not
         % matter theoretically.        
        %g = grad(prob,x_cg)
        %grad_fullspace = grad_fullspace - (x_cg.U*g.M*x_cg.V' + x_cg.U*g.Vp' + g.Up*x_cg.V');
        
        if opts.rand_pca
            [U,S,V] = pca(grad_fullspace/p,rank_increase,true,4);%,10,2*rank_increase+5); 
        else
            [U,S,V,flag] = svds(grad_fullspace/p, rank_increase, 'L', opts_svds);        
            if flag, warning('svds did not converge'); end        
        end
        
        
        %U = U - x_cg.U*(x_cg.U'*U);
        %V = V - x_cg.V*(x_cg.V'*V);
        %subspace(U, x_cg.U)
        %subspace(V, x_cg.V)
        
        % do line search and  make new x
        %tmin = exact_search_onlyTxM_LR(prob,x_cg,U*S,V);
        
        %--> try tmin = -1 without opt.S, doe opt.S lals optie
        %tmin = -1;  % also works
    else
        U = orth(rand(n, rank_increase));
        V = orth(rand(n, rank_increase));
        U = U - x_cg.U*(x_cg.U'*U);
        V = V - x_cg.V*(x_cg.V'*V);
        S = diag(rand(rank_increase,1)); S = x_cg.sigma(end)*S/norm(S,2);                
    end
    
    x0.U = [x_cg.U U];
    x0.V = [x_cg.V V];
    
    if opts.opt_S    
        x0.U = orth(x0.U);
        x0.V = orth(x0.V);        
        
        %tmin = exact_search_onlyTxM_LR(prob,x_cg,U*S,V);
        %S0 = blkdiag(diag(x_cg.sigma), tmin*S);
        
%         TT = tic();
         S = getoptS_2(prob,x0.U,x0.V,M_omega,opts);                     
%         oS = toc(TT)
%        TT = tic();
%        S = getoptS_4(prob,x0.U,x0.V,M_omega_old,opts);                     
%        oS_2 = toc(TT)
        
        
    else
        tmin = exact_search_onlyTxM_LR(prob,x_cg,U*S,V);
        %if opts.increase_with_residual
        %    tmin = -1;
        %end
        S = blkdiag(x_cg.S, tmin*S);
    end
    
    
    
    [x0.U, Ru] = qr(x0.U, 0);
    [x0.V, Rv] = qr(x0.V, 0);
    [UU,SS,VV] = svd(Ru*S*Rv');
    x0.S = SS;
    x0.U = x0.U*UU;
    x0.V = x0.V*VV;

    ranklist = [ranklist,prob.r];
    outs_c.N = N_c;
    if isfield(opts_cg,'saveiterates') && opts_cg.saveiterates == 1
        Xr_c = cell(1,outs_c.N);
        for kk=1:outs_c.N
            outs_c.Xout{kk} = outs_c.Xout{kk}';
            Xr_c{kk}={outs_c.Xout{kk}{1},outs_c.Xout{kk}{2}}';
        end
        outs_c.Xout = outs_c.Xout';
    else
        outs_c.Xout{1}=outs_c.Xout{1}';
        Xr_c=outs_c.Xout;
    end
    outs_c.U={outs_c.U};
    outs_c.V={outs_c.V};
    outs_c.S={outs_c.S};
    outs_c.on_omega={outs_c.on_omega};
    if prob.r == rank_increase
        outs_c.time = tm_c';
        outs = outs_c;
        if isfield(opts_cg,'saveiterates') && opts_cg.saveiterates == 1
            Xr   = Xr_c';
        else
            Xr   = Xr_c;
        end
    else
        outs_c.time = tm_c'+outs.time(end);
        outs = cell2struct(cellfun(@vertcat,struct2cell(outs),...
            struct2cell(outs_c),'uni',0),fieldnames(outs),1);
        if isfield(opts_cg,'saveiterates') && opts_cg.saveiterates == 1
            Xr = [Xr;Xr_c'];
        else
            Xr = [Xr;Xr_c];
        end
        %Xr = cell2struct(cellfun(@vertcat,struct2cell(Xr),...
          %  struct2cell(Xr_c),'uni',0),fieldnames(Xr),1);
    end
    prob.r = prob.r + rank_increase;
    x0 = prepx(prob, x0);
    its = its+1;
end
fprintf('Total number of inner iterations = %4i, total time = %5.2e (sec) \n', size(hist,1), toc(t))
Xr=Xr';
outs.ranklist = ranklist;
outs.N_vec = outs.N;
outs.N = sum(outs.N_vec);
x = x0;
histout = hist;
stats.time = toc(t);
stats.chg_ranks = chg_ranks;

end

function S = getoptS(prob,L,R, M_omega,opts)

% Solve projected LSQ via QR. Turns out to be much slower than via normal
% equations.

[n, r] = size(L);
[m, r] = size(R);

% C = L' * (M_omega*R);
% C = C(:) ; %right-hand side vector

if opts.opt_S==2
    C = prob.opt_S_fast.data;
else
    C = prob.data;
end



A = zeros(length(C), r^2);
for i = 1:r
    for j = 1:r
        ind = (j-1)*r + i ;
        if opts.opt_S==2
            d = partXY(L(:,i)', R(:,j)', prob.opt_S_fast.Omega_i, prob.opt_S_fast.Omega_j, prob.opt_S_fast.m)';
            %updateSval(prob.opt_S_fast.temp_omega, d, length(d));
            %temp = L' * (prob.opt_S_fast.temp_omega * R);
        elseif opts.opt_S==1
            d = partXY(L(:,i)', R(:,j)', prob.Omega_i, prob.Omega_j, prob.m)';
            %updateSval(prob.temp_omega, d, length(d));
            %temp = L' * (prob.temp_omega * R);
        else
            error('opt_S')
        end
        temp = d;
        A(:,ind) = temp(:);
    end
end

S = A\C ;
S = reshape(S,r,r) ;
end

function S = getoptS_2(prob,L,R, M_omega,opts)

[n, r] = size(L);
[m, r] = size(R);

C = L' * (M_omega*R);
C = C(:) ; %right-hand side vector

A = zeros(r^2);
for i = 1:r
    for j = 1:r
        ind = (j-1)*r + i ;
        if opts.opt_S==2
            d = partXY(L(:,i)', R(:,j)', prob.opt_S_fast.Omega_i, prob.opt_S_fast.Omega_j, prob.opt_S_fast.m)';
            updateSval(prob.opt_S_fast.temp_omega, d, length(d));
            temp = L' * (prob.opt_S_fast.temp_omega * R);
        elseif opts.opt_S==1
            d = partXY(L(:,i)', R(:,j)', prob.Omega_i, prob.Omega_j, prob.m)';
            updateSval(prob.temp_omega, d, length(d));
            temp = L' * (prob.temp_omega * R);
        else
            error('opt_S')
        end        
        A(:,ind) = temp(:);
    end
end

S = A\C ;
S = reshape(S,r,r) ;
end

function S = getoptS_3(prob,L,R, opts)
% does not work: the subsampled normal equations are not a good
% preconditioner!


[n, r] = size(L);
[m, r] = size(R);

A = zeros(prob.opt_S_fast.m,r^2);
for i = 1:r
    for j = 1:r
        ind = (j-1)*r + i ;
        
        d = partXY(L(:,i)', R(:,j)', prob.opt_S_fast.Omega_i, prob.opt_S_fast.Omega_j, prob.opt_S_fast.m)';
        %updateSval(prob.opt_S_fast.temp_omega, d, length(d));
        %temp = L' * (prob.opt_S_fast.temp_omega * R);
        
        A(:,ind) = d(:);
    end
end

[Q,prec] = qr(A,0);


    function y = fun_A(x,transp_flag)
        if strcmp(transp_flag,'transp')
            updateSval(prob.temp_omega, x, length(x));
            Y = L'*(prob.temp_omega * R);
            y = Y(:);
          
       elseif strcmp(transp_flag,'notransp')
          S = reshape(x,r,r);
          y = partXY(S'*L', R', prob.Omega_i, prob.Omega_j, prob.m)';
       end
    end

%prec = [];
tol = 1e-5;
[x,flag,res,iter] = lsqr(@fun_A,prob.data,tol,100,prec);

iter
S = reshape(x,r,r) ;
end

function S = getoptS_4(prob,L,R, M_omega, opts)
% does not work: the subsampled normal equations are not a good
% preconditioner!
[n, r] = size(L);
[m, r] = size(R);
C = L' * (M_omega*R);
C = C(:);

A = zeros(r^2);
for i = 1:r
    for j = 1:r
        ind = (j-1)*r + i ;
        
        d = partXY(L(:,i)', R(:,j)', prob.opt_S_fast.Omega_i, prob.opt_S_fast.Omega_j, prob.opt_S_fast.m)';
        updateSval(prob.opt_S_fast.temp_omega, d, length(d));
        temp = L' * (prob.opt_S_fast.temp_omega * R);
        
        A(:,ind) = temp(:);
    end
end


    function y = fun_A(x)     
        S = reshape(x,r,r);
        y = partXY(S'*L', R', prob.Omega_i, prob.Omega_j, prob.m)';        
        updateSval(prob.temp_omega, y, length(y));
        temp = L' * (prob.temp_omega * R);
        y = temp(:);
    end

prec = A;
tol = 1e-5;
[x,flag,res,iter] = pcg(@fun_A,C,tol,100,prec);

iter
S = reshape(x,r,r) ;
end
