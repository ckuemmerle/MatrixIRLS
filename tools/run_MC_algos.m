function [Xr,outs,alg_name] = run_MC_algos(Phi,y,r,alg_name,opts_custom)
% This function runs different algorithms (indicated by alg_name)
% for a given matrix completion problem with entry mask Phi and provided 
% data y.

% =========================================================================
% Parameters
% ----------
% Phi:  (d1 x d2) sparse matrix. Non-zero indices correspond to the mask of
%       known matrix entries.
% y:    (m x 1) vector of samples. Provided matrix entry values,
%       ordered in according to linear indices of non-zero entries of Phi.
% r:    Target rank to be used.
% alg_name:     (1 x nr_algos) cell of character strings. Indicating which 
%               algorithm to use.
% opts_custom:  struct. Options to override standard options of any 
%               algorithm.
% Returns
% ----------
% Xr:   (1 x nr_algos_total) cell. The l-th cell contains algorithmic
%       iterates of l-th algorithm (if opts.saveiterates == 1) or
%       the last iterate (if opts.saveiterates == 0).
% outs: (1 x nr_algos_total) cell. The l-th cell contains additional
%       information about the progress of the l-th algorithm.
% alg_name: (1 x nr_algos_total) cell of character strings. Contains name
%       of used algorithms, including reflecting potentially additional 
%       parameter choices (and thus, might differ from input 'alg_name').
% =========================================================================
% Author: Christian Kuemmerle, 2018-2020.

[d1,d2]=size(Phi);
Omega = find(Phi);
nr_algos = length(alg_name);
opts_new = cell(1,nr_algos);
if isfield(opts_custom,'r')
   r = opts_custom.r;
   opts_custom=rmfield(opts_custom,'r');
end

%%% prepare multiple algorithm instances of IRLS variants if multiple
%%% values for quasi-norm parameter p are provided
if isfield(opts_custom,'p')
    len_p    = length(opts_custom.p);
    if len_p >= 1
        ind=find(contains(alg_name,'HM-IRLS')+contains(alg_name,'AM-IRLS')+...
            contains(alg_name,'IRLS-col')+contains(alg_name,'IRLS-row')+...
            contains(alg_name,'HM_IRLS_PCG')+contains(alg_name,'HM_IRLS_concave')+...
            contains(alg_name,'MatrixIRLS')+contains(alg_name,'IRLS-M')+...
            contains(alg_name,'sIRLS-p')+contains(alg_name,'IRLS-p')+...
            contains(alg_name,'IRucLq')+contains(alg_name,'tIRucLq'));
        for i=1:length(ind)
            ind_c=ind(i)+(i-1)*(len_p-1);
            alg_name_ind_c=cell(1,len_p);
            for ii=1:len_p
                alg_name_ind_c{ii}=[alg_name{ind_c},' p=',num2str(opts_custom.p(ii))];
            end
            if ind_c < length(alg_name)
                alg_name= [alg_name(1:ind_c-1),alg_name_ind_c,alg_name(ind_c+1:end)];
                opts_new= [opts_new(1:ind_c-1),repelem(opts_new(ind_c),len_p),opts_new(ind_c+1:end)];
            else
                alg_name= [alg_name(1:ind_c-1),alg_name_ind_c];
                opts_new= [opts_new(1:ind_c-1),repelem(opts_new(ind_c),len_p)];
            end
            for j=1:len_p
                opts_new{ind_c+j-1}.p=opts_custom.p(j);
            end
        end
    end
end

if isfield(opts_custom,'type_mean')
    len_typmean = length(opts_custom.type_mean);
    if len_typmean >= 1
        ind=find(contains(alg_name,'MatrixIRLS'));
        for i=1:length(ind)
            ind_c=ind(i)+(i-1)*(len_typmean-1);
            alg_name_ind_c=cell(1,len_typmean);
            for ii=1:len_typmean
                alg_name_ind_c{ii}=[alg_name{ind_c},' meantype=',num2str(opts_custom.type_mean{ii})];
            end
            if ind_c < length(alg_name)
                alg_name= [alg_name(1:ind_c-1),alg_name_ind_c,alg_name(ind_c+1:end)];
                opts_new= [opts_new(1:ind_c-1),repelem(opts_new(ind_c),len_typmean),opts_new(ind_c+1:end)];
            else
                alg_name= [alg_name(1:ind_c-1),alg_name_ind_c];
                opts_new= [opts_new(1:ind_c-1),repelem(opts_new(ind_c),len_typmean)];
            end
            for j=1:len_typmean
                opts_new{ind_c+j-1}.type_mean=opts_custom.type_mean{j};
            end
        end
    end
end

nr_algos = length(alg_name);
names = fieldnames(opts_custom);
for l=1:nr_algos
    for i=1:length(names)
        if not(isequal(names(i),{'p'}))
            if not(iscell(opts_custom.(names{i})))
                opts_new{l}.(names{i})=opts_custom.(names{i});
            end
        end
    end
end



prob  = cell(1,nr_algos);
outs  = cell(1,nr_algos);
Xr    = cell(1,nr_algos);
opts  = cell(1,nr_algos);
lambda = cell(1,nr_algos);
params = cell(1,nr_algos);

%%% run algorithm (with default options, priority on option in opts_new if
%%% applicable)
for l=1:nr_algos
    if contains(alg_name{l},'HM-IRLS') || contains(alg_name{l},'AM-IRLS') ...
            || contains(alg_name{l},'IRLS-col') || contains(alg_name{l},'IRLS-row')
        opts{l} = getDefaultOpts_IRLS;
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l}.linsolve = 'direct_solve';
        prob{l}.d1     = d1;
        prob{l}.d2     = d2;
        prob{l}.r      = r;
        prob{l}.Omega  = Omega;
        prob{l}.Phi    = Phi;
        prob{l}.y      = y;
        if contains(alg_name{l},'HM-IRLS')
                [Xr{l},N,eps,singval,time] = ...
            HM_IRLS(prob{l},opts{l});
        elseif contains(alg_name{l},'AM-IRLS')
                [Xr{l},N,eps,singval,time] = ...
            AM_IRLS(prob{l},opts{l});
        elseif contains(alg_name{l},'IRLS-col')
            opts{l}.variant='left';
                [Xr{l},N,eps,singval,time] = ...
            IRLS_onesided(prob{l},opts{l});
        elseif contains(alg_name{l},'IRLS-row')
            opts{l}.variant='right';
                [Xr{l},N,eps,singval,time] = ...
            IRLS_onesided(prob{l},opts{l});  
        end
        outs{l}.N   = N;
        outs{l}.eps = eps;
        outs{l}.singval = singval;
        outs{l}.time    = time;
    elseif contains(alg_name{l},'MatrixIRLS')        
        opts{l} = getDefaultOpts_IRLS;
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        prob{l}.d1     = d1;
        prob{l}.d2     = d2;
        prob{l}.r      = r;
        prob{l}.Phi    = Phi;
        if iscell(y) && opts{l}.recsys
            prob{l}.y       =y{1};
            prob{l}.test    =y{2};
            prob{l}.train   =y{3};
            prob{l}.centscale = y{4};
        else
            prob{l}.y      = y;
        end
        if isfield(opts{l},'lambda') && ~isempty(opts{l}.lambda)
            lambda{l} = opts{l}.lambda;
        else
            lambda{l} = 0;
        end
        [X_c,outs{l}] = ...
        MatrixIRLS(prob{l},lambda{l},opts{l});
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            for kk=1:outs{l}.N 
                Xr{l}{kk}=outs{l}.X{kk};
            end
        else
             Xr{l}{1} = X_c;
        end
    elseif contains(alg_name{l},'sIRLS-p') || contains(alg_name{l},'IRLS-p')
        %%% Algorithms by [Mohan & Fazel (2012)]
        opts{l} = getDefaultOpts_IRLS;
        prob{l}.d1     = d1;
        prob{l}.d2     = d2;
        prob{l}.r      = r;
        prob{l}.Phi    = Phi;
        [rowind,colind]=find(Phi);
        prob{l}.X0_revealed=sparse(rowind,colind,y,d1,d2);
        
        opts{l}.gam0 = 1e-2;
        opts{l}.gammin = 1e-10;
        if not(isempty(r))
            fr = r*(d1+d2-r)/length(y);
            if(fr < 0.4)
                opts{l}.eta = 1.1;
                opts{l}.incr = 50;
            else
                opts{l}.eta = 1.03;
                opts{l}.incr = 100;
            end
        else
        	opts{l}.eta = 1.03;
        	opts{l}.incr = 100;
        end
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        if contains(alg_name{l},'sIRLS-p')
            [Xr_c,outs{l}] = sirls_p_adp(prob{l},opts{l});
        else
            [Xr_c,outs{l}] = irls_p_adp(prob{l},opts{l});
        end
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            Xr{l} = cell(1,outs{l}.N);
            for kk=1:outs{l}.N
                Xr{l}{kk}=outs{l}.X{kk};
            end
        else
            Xr{l} = {Xr_c};
        end
    elseif contains(alg_name{l},'IRucLq') || contains(alg_name{l},'tIRucLq')
        %%% Algorithms by [M.-J. Lai, Y. Xu, W. Yin (2013)]
        opts{l} = getDefaultOpts_IRLS;
        prob{l}.d1     = d1;
        prob{l}.d2     = d2;
        prob{l}.r      = r;
        prob{l}.Phi    = Phi;
        [rowind,colind]=find(Phi);
        prob{l}.X0_revealed=sparse(rowind,colind,y,d1,d2);
        opts{l}.gamma = 0.9;
        opts{l}.gamma0 = 0.9;
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        if opts{l}.lambda == 0
            lambda_c = 1e-6.*norm(y);
        else
            lambda_c = opts{l}.lambda;
        end
        if contains(alg_name{l},'tIRucLq')
            [Xr_c,outs{l}] = tIRucLq_m_adp(prob{l},lambda_c,opts{l});
        else
            [Xr_c,outs{l}] = IRucLq_m_adp(prob{l},lambda_c,opts{l});
        end
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            Xr{l} = cell(1,outs{l}.N);
            for kk=1:outs{l}.N
                Xr{l}{kk}=outs{l}.X{kk};
            end
        else
            Xr{l} = {Xr_c};
        end
    elseif strcmp(alg_name{l},'LMaFit')
        opts{l}.tol=1.25e-10;
        opts{l}.print=0;
        opts{l}.est_rank=0; 
        opts{l}.Zfull=1;
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l}.maxit=opts{l}.N0_firstorder; 
        if iscell(y) && opts{l}.recsys
            opts{l}.test    =y{2};
            opts{l}.train   =y{3};
            opts{l}.centscale = y{4};
            y_LMa           =y{1};
        else
            y_LMa           =y;
        end
        [~,~,Out] = lmafit_mc_adp(d1,d2,r,Omega,y_LMa,opts{l});
        outs{l}=Out;
        outs{l}.N=Out.iter;
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            for kk=1:outs{l}.N
                Xr{l}{kk}={Out.Xhist{kk},Out.Yhist{kk}'};
            end
        else
            Xr{l}{1}={Out.Xhist,Out.Yhist'};
        end
       
   elseif strcmp(alg_name{l},'LRGeomCG')
       opts{l} = setExtraOpts(opts{l},opts_new{l});
        if iscell(y) && opts{l}.recsys
            prob{l} = make_prob_Riem_adp(y{1},Omega,d1,d2,r);
            prob{l}.test    =y{2};
            prob{l}.train   =y{3};
        else
            prob{l} = make_prob_Riem_adp(y,Omega,d1,d2,r);
        end
        start = make_start_x_Riem(prob{l});
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l} = default_opts_Riem(opts{l}.N0_firstorder,opts{l}.tol);
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l}.N0 = opts{l}.N0_firstorder;
        [~,~,~,outs{l},N,tm_c] = LRGeomCG_outp(prob{l},opts{l},start);
        outs{l}.N=N;
        outs{l}.time=tm_c;
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            for kk=1:outs{l}.N
                Xr{l}{kk}={outs{l}.Xout{kk}{1},outs{l}.Xout{kk}{2}};
            end
        else
            Xr{l}=outs{l}.Xout;
        end
    elseif strcmp(alg_name{l},'ScaledGD')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        prob{l}.d1 = d1;
        prob{l}.d2 = d2;
        opts{l}.T = opts{l}.N0_firstorder;
        opts{l}.eta = 0.5; 
        [prob{l}.rowind,prob{l}.colind]=find(Phi);
        [prob{l}.Omega] = find(Phi);
        prob{l}.y = y;
        [Xr{l},outs{l}] = ScaledGD(prob{l},opts{l},r);
    elseif strcmp(alg_name{l},'BFGD')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        [rowind,colind]=find(Phi);
        if iscell(y) && opts{l}.recsys
            prob{l}.y         =y{1};
            [f_c,f_grad_c,params{l}] = BFGD_aux(d1,d2,r,...
          rowind,colind,Omega,prob{l}.y,1,opts{l}.N0_firstorder);
            params{l}.test    =y{2};
            params{l}.train   =y{3};
            params{l}.centscale = y{4};
            params{l}.clip_ratingscale = opts{l}.clip_ratingscale;
        else
            prob{l} = make_prob_Riem_adp(y,Omega,d1,d2,r);
            prob{l}.y                  = y;
            [f{l},f_grad{l},params{l}] = BFGD_aux(d1,d2,r,...
            rowind,colind,Omega,prob{l}.y,1,opts{l}.N0);
        end
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            params{l}.saveiterates = 1;
        end
        if isfield(opts{l},'recsys') && opts{l}.recsys == 1
            params{l}.recsys = 1;
        end
        X0_revealed=sparse(rowind,colind,prob{l}.y,d1,d2);
        tic
        [Xr{l},outs{l}] = BFGD_outp(f_c,f_grad_c,params{l},X0_revealed);
        outs{l}.N=outs{l}.numiter;
        outs{l}.time=toc;
    elseif strcmp(alg_name{l},'NNM')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        outs{l} = struct;
        outs{l}.N=1;
        Phi_sorted = sparse(1:length(Omega),Omega,ones(length(Omega),1),length(Omega),d1*d2);
        tic
        TEST=nuclear_norm_min(d1,d2,Phi_sorted,y,1);
        Xr{l} = {TEST};
        alg_name{l}='NNM';
        outs{l}.time=toc;
    elseif strcmp(alg_name{l},'NNM_DouglasRachford')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        tic
%         [X{i}(:,:,1)] = nuclear_norm_min(d1,d2,Phi_sorted,aux_indic(:,3),1);
        outs{l} = struct;
        [Xr{l},outs{l}.N,~] = NNM_DouglasRachford(d1,d2,Omega,y,opts{l}.N0);
        alg_name{l}='NNM';
        outs{l}.time=toc;
    elseif strcmp(alg_name{l},'NIHT')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        start=make_start_x_IHT_ASD(alg_name{l},d1,d2,r,Omega,y);
%         opts{i} = default_opts();
        opts{l}= opts_adapted(opts{l}.N0_firstorder,opts{l}.tol,opts{l}.tol,...
            opts{l}.verbose,opts{l}.saveiterates);
        [Xr{l},outs{l}] = NIHT_Matrix_outp(d1,d2,r,Omega,y,start,opts{l});
        outs{l}.N=outs{l}.iter; 
    elseif strcmp(alg_name{l},'ASD')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        start = make_start_x_IHT_ASD(alg_name{l},d1,d2,r,Omega,y);
%         opts{l} = default_opts();
        opts{l} = opts_adapted(opts{l}.N0_firstorder,opts{l}.tol,opts{l}.tol,...
            opts{l}.verbose,opts{l}.saveiterates);
        [Xr{l},outs{l}] = ASD_outp(d1,d2,r,Omega,y,start,opts{l});
        outs{l}.N=outs{l}.iter;
    elseif strcmp(alg_name{l},'ScaledASD')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        start = make_start_x_IHT_ASD(alg_name{l},d1,d2,r,Omega,y);
%         opts{l} = default_opts();
        opts{l}= opts_adapted(opts{l}.N0_firstorder,opts{l}.tol,opts{l}.tol,...
            opts{l}.verbose,opts{l}.saveiterates);
        [Xr{l},outs{l}] = ScaledASD_outp(d1,d2,r,Omega,y,start,opts{l});
        outs{l}.N=outs{l}.iter; 
    elseif strcmp(alg_name{l},'OptSpace')
        [rowind,colind]=find(Phi);
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        X0_revealed=sparse(rowind,colind,y,d1,d2);
        [~,~,~,~,Xr{l},outs{l}.N,outs{l}.time] = OptSpace_output(X0_revealed,...
            r,opts{l}.N0_firstorder,opts{l}.tol,opts{l}.verbose);     
    elseif strcmp(alg_name{l},'R2RILS') || strcmp(alg_name{l},'R2RILS_svd')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l}.assemblesparseA = 1;
        opts{l}.N0 = 30;
        matind = zeros(length(Omega),2);
        [matind(:,1), matind(:,2)] = ind2sub([d1,d2],Omega); 
        Y_mat=sparse(matind(:,1),matind(:,2),y,d1,d2);
        if strcmp(alg_name{l},'R2RILS')
            [Xr{l},outs{l}] = R2RILS_adp(Y_mat,matind,r,opts{l});
        else
            [Xr{l},outs{l}] = R2RILS_svd(Y_mat,matind,r,opts{l});
        end
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            Xr{l} = cell(1,outs{l}.N);
            for kk=1:outs{l}.N
                Xr{l}{kk}=outs{l}.Xout{kk};
            end
        else
            Xr_c = Xr{l};
            Xr{l} = {Xr_c};
        end
    elseif strcmp(alg_name{l},'RTRMC')
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l}.method = 'rtr';     % 'cg' or 'rtr', to choose the optimization algorithm
        opts{l}.order = 2;          % for rtr only: 2 if Hessian can be used, 1 otherwise
        opts{l}.precon = true;      % with or without preconditioner
        opts{l}.maxiter = opts{l}.N0;      % stopping criterion on the number of iterations
        opts{l}.storedepth = opts{l}.maxiter;
        opts{l}.maxinner = opts{l}.N0_inner;      % for rtr only : maximum number of inner iterations
        opts{l}.tolgradnorm = opts{l}.tol*1e-6;%1e-8; % stopping criterion on the norm of the gradient
        opts{l}.verbosity = opts{l}.verbose;      % how much information to display during iterations
        opts{l}.computeRMSE = false; % set to true if RMSE is to be computed at each step
        if isfield(opts{l},'lambda') && ~isempty(opts{l}.lambda)
            lambda{l}=opts{l}.lambda;
        else
            lambda{l}= 0;
        end
        lambda{l}=max(1e-8,lambda{l});
        [Om_col,Om_row] = ind2sub([d1,d2],Omega); 
        C = ones(length(Omega),1); % confidence vector
        prob{l} = buildproblem(Om_col, Om_row, y, C, d1, d2, r, lambda{l});
        initstart = tic;
        U0 = initialguess(prob{l});
        inittime = toc(initstart);
        [Uc, Wc, stats_c] = rtrmc(prob{l}, opts{l}, U0);
        for i = 1 : length(stats_c)
            stats_c(i).time = stats_c(i).time + inittime;
            if i > 1
                outs{l}.time(i-1) = stats_c(i).time;
            end
        end
        
        iter_tmp = [stats_c.iter];
        N_c = iter_tmp(end);
        outs{l}.N = N_c;
        outs{l}.stats = stats_c;
        if isfield(opts{l},'saveiterates') && opts{l}.saveiterates == 1
            Xr{l} = cell(1,outs{l}.N);
            for kk=1:outs{l}.N
                Xr{l}{kk}={stats_c(kk+1).U,stats_c(kk+1).W.'};
            end
        else
            Xr{l}={{Uc,Wc.'}};
        end
    end
    outs{l}.opts=opts{l};
end


end

function opts = setExtraOpts(opts,opts_new)
    if ~isempty(opts_new)
        optNames = fieldnames(opts_new);
        for i = 1:length(optNames)
            optName = optNames{i};
            opts.(optName) = opts_new.(optName);
        end
    end

end

