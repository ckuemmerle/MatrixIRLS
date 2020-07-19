function [X,Y,Out] = lmafit_mc_adp(m,n,k,Known,data,opts)
%
% Implementation of "Low-rank Matrix Fitting algorithm" ('LMaFit') for
% matrix completion [1].
% Code by Yin Zhang and Zaiwen Wen.
% =========================================================================
% [1] Zaiwen Wen, Wotao Yin, and Yin Zhang. Solving a low-rank 
% factorization model for matrix completion by a nonlinear successive 
% overrelaxation algorithm, Mathematical Programming Computation, vol. 4, 
% no. 4, pp. 333?361, 2012.
% =========================================================================
% Solver for matrix completion:
%
%        U_ij ~= A_ij,  for i,j in Known
%
% where U =  X*Y and A_ij's are given in an known index set.
%
% Output:
%           X --- m x k matrix
%           Y --- k x n matrix
%         Out --- output information
% Input:
%        m, n --- matrix sizes
%           k --- rank estimate
%        data --- values of known elements in a 1D row vector
%       Known --- positions of known elements in a 1D row vector
%                 assuming matrices are arranged column-wise
%             or  Known is a structure with two fields [Ik, Jk]
%        opts --- option structure with fields as follows:
%
%        ******************************************************************
%        Fields in opts:
%        tol      -    tolerance, default: 1.25e-4
%        maxit    -    max. number of iterations, default:  500
%        print    = {0, 1, 2}; level of output information
%        DoQR     = 1; use QR to solve the least square subproblem
%                 = 0; use "linsolve" ....
%        Zfull    = 0; use a dense matrix to store Z
%                   1; use the factor and sparse form Z = XY + S
%        est_rank = 0; rank is fixed
%                 = 1; the decreasing rank strategy, 
%                      related fields: rank_min, rk_jump
%                 = 2; the increasing rank strategy
%                      related fields: rank_max, rk_inc
%        rank_min -    the minimal rank allowed if est_rank = 1
%                      default: 1
%        rank_max -    the maximal rank allowed if est_rank = 2
%                      default: max(floor(0.1*min(m,n)),2*k)
%        rk_jump  -    tolerance to decrease the rank if est_rank = 1
%                      default: 10
%        rk_inc   -    the increase of the rank if est_rank = 2
%                      default: 1; 
%					   *Other value can be max(min(5,floor(rank_max/5)),1);
%					   *rk_inc is doubled if the current rank estimate rk
%                       is larger than 50
%        init     = 0; no initial solutions is provided
%                 = 1; initial solutions are stored in opt.X and opts.Y; 
%        ******************************************************************
%
% Copyright(c) 2009 Yin Zhang
%
%   modified by Zaiwen Wen, 12/17/2009
%   modified by Zaiwen Wen, 09/30/2011

L = length(data);

% set parameters
tol = 1.25e-4;
maxit = 500;
iprint = 1;
Zfull = (L/(m*n) > 0.2 ) || k > .02*min(m,n) || m*n < 5e5;
DoQR = true;
est_rank = 1;
rank_max =  max(floor(0.1*min(m,n)),2*k);
rank_min =  1;
rk_inc = 1; %max(min(5, floor(rank_max/5)),1);
rk_jump = 10;
init = 0;
save_res = 0;

if isfield(opts,'tol');         tol     = opts.tol;        end
if isfield(opts,'maxit');       maxit   = opts.maxit;      end
if isfield(opts,'print');       iprint  = opts.verbose;      end
if isfield(opts,'Zfull');       Zfull   = opts.Zfull;      end
if isfield(opts,'DoQR');        DoQR    = opts.DoQR;       end
if isfield(opts,'est_rank');    est_rank= opts.est_rank;   end
if isfield(opts,'rank_max');    rank_max= opts.rank_max;   end
if isfield(opts,'rank_min');    rank_min= opts.rank_min;   end
if isfield(opts,'rk_inc');      rk_inc  = opts.rk_inc;     end 
if isfield(opts,'rk_jump');     rk_jump = opts.rk_jump;    end 
if isfield(opts,'init');        init    = opts.init;       end
if isfield(opts,'save_res');    save_res= opts.save_res;   end

if isfield(opts,'recsys') && opts.recsys == 1
    recsys = 1;
    if isfield(opts,'train') && isfield(opts,'test')
    else
        error('Please provide training and test data (recommender system application).')
    end
    clip_ratingscale = opts.clip_ratingscale;
    if clip_ratingscale
        clip_ratingscale = [min(opts.train.values), max(opts.train.values)];
    end
    RMSE = zeros(1,maxit);
    RMSE_train = zeros(1,maxit);
    MABS = zeros(1,maxit);
    MABS_train = zeros(1,maxit);
else
    recsys = 0;
end
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    Out.Xhist=cell(1,maxit);
    Out.Yhist=cell(1,maxit);
end

reschg_tol = 0.5*tol; rk = k;  
if est_rank == 1; rank_max = min(rank_max, k); end 
    
linopts.SYM = true; linopts.POSDEF = true;
datanrm = max(1,norm(data));    

objv = zeros(maxit,1);  RR = ones(maxit,1);
if iprint == 1; fprintf('Iteration:     '); end
if iprint == 2
    fprintf('\nLMafit_mc: Zfull = %i, DoQR = %i\n',Zfull,DoQR);
end

% initialize: make sure the correctness of the index set and data
data(data==0) = eps;    data_tran = false; 


if Zfull %Z is full
    if isstruct(Known), Known = sub2ind([m n], Known.Ik,Known.Jk); end
     Z = zeros(m,n);
%     size(Known)
%     size(data)
    Z(Known) = data';   
else %Z = S + XY, initialize the storage of S
    if isnumeric(Known);       [Ik,Jk] = ind2sub([m n],Known);
    elseif isstruct(Known)     Ik = Known.Ik; Jk = Known.Jk;     end
    %make sure the order of Ik, Jk and data are correctly as in a sparse matrix
    S = sparse(Ik, Jk, data, m, n);  [Ik, Jk, data] = find(S);  data = data';   
end

if m>n 
    tmp = m; m = n; n = tmp; data_tran = true; 
    if Zfull;  Z = Z';  Known = find(Z);        data = Z(Known);
    else       S = S';  [Ik, Jk, data] = find(S);  data = data';    end
    if recsys 
        tmp=opts.test.col;
        opts.test.col=opts.test.row;
        opts.test.row=tmp;
        tmp=opts.train.col;
        opts.train.col=opts.train.row;
        opts.train.row=tmp;
    end
end

if init == 0
    X = zeros(m,k);   Y = eye(k,n);   Res = data;   res = datanrm; 
    %X = zeros(m,k);   Y = rand(k,n);   Res = data;   res = datanrm; 
else
    X = opts.X;  Y = opts.Y;    opts.X = []; opts.Y = [];
    if Zfull 
        Z = X*Y;  Res = data - Z(Known);    Z(Known) = data;
    else %Z = S + XY, initialize the storage of S
        Res = data - partXY(X',Y,Ik,Jk,L);  updateSval(S, Res, L);
    end
    res = norm(Res);
end

% parameters for alf
alf = 0;  increment = 1; %init_inc();
itr_rank = 0; minitr_reduce_rank = 5;    maxitr_reduce_rank = 50;  

time=zeros(1,maxit);
tic
% main iteration
for iter = 1:maxit
    itr_rank = itr_rank + 1; 
    Xo = X; Yo = Y; Res0 = Res; res0 = res; alf0x = alf;
    if Zfull
        Zo = Z; X = Z*Y';
        if est_rank == 1
            [X,R,E] = qr(X,0);  Y = X'*Z;
        elseif DoQR
            [X,R  ] = qr(X,0);  Y = X'*Z;
        else
            Xt = X'; Y = linsolve(Xt*X,Xt*Z,linopts);
        end
%         size(Res)
%         size(Z(Known))
        Z = X*Y;     Res = data - Z(Known); 
    else % Z=S+XY for sparse S(Known)=data-XY(Known)
        Yt = Y';     X = S*Yt + X*(Y*Yt);    % Z*Y'
        if est_rank == 1
            [X,R,E] = qr(X,0);  Xt = X';
            Y = Xt*S + (Xt*Xo)*Y;       % X'*Z
        elseif DoQR
            [X,R  ] = qr(X,0);  Xt = X';
            Y = Xt*S + (Xt*Xo)*Y;
        else
            Xt = X';
            Y = Xt*S + (Xt*Xo)*Y;
            Y = linsolve(Xt*X,Y,linopts);
        end
        Res = data - partXY(Xt,Y,Ik,Jk,L);
    end %Zfull

    res = norm(Res);  relres = res/datanrm;  ratio = res/res0;
    reschg = abs(1-res/res0); RR(iter) = ratio; 
    
    % rank estimation. if rank is reduced, res can be increased
    if est_rank >= 1; rank_estimator_adaptive(); end
    if rk ~= k; k = rk; if est_rank ==0; alf = 0; continue; end; end

    % adjust alf
    if ratio >= 1 
        increment = max(0.1*alf, 0.1*increment);
        %increment = max(0.1*alf, 0.5*increment);
        X = Xo; Y = Yo; Res = Res0; res = res0; relres = res/datanrm;
        alf = 0;    if Zfull; Z = Zo;end
    elseif ratio > 0.7
        increment = max(increment, 0.25*alf);   
        alf = alf + increment;
    end

    % printout
    if iprint == 1; fprintf('\b\b\b\b\b%5i',iter); end
    if iprint == 2
        fprintf('it: %5i rk: %4d, rel. %3.1e r. %4.4f chg: %3.1e alf: %3.1e inc: %3.1e\n',...
            iter, k, relres,ratio,reschg,alf0x,increment);
    end
    objv(iter) = relres;
    time(iter)=toc;

    % update Z or S with the most recent alf
    if Zfull; Z(Known) = data + alf*Res; else updateSval(S, (alf+1)*Res, L); end
    if isfield(opts,'saveiterates') && opts.saveiterates == 1
        if m>n
            Out.Xhist{iter}=X;
            Out.Yhist{iter}=Y;
        else
            Out.Xhist{iter}=X;%Y';
            Out.Yhist{iter}=Y;%X';
        end
    end
    if recsys 
        Yt=Y';
        XX=X;
        prediction_train_val_c =  get_prediction_RecSys_MatFac(XX, Yt, ...
            eye(size(XX,2),size(Yt,2)),opts.train.rowind, opts.train.colind,...
            opts.centscale,clip_ratingscale);
        prediction_test_val_c =  get_prediction_RecSys_MatFac(XX, Yt, ...
            eye(size(XX,2),size(Yt,2)),opts.test.rowind, opts.test.colind,...
            opts.centscale,clip_ratingscale);
        RMSE(iter) = get_RMSE(prediction_test_val_c,opts.test.values);
        RMSE_train(iter) = get_RMSE(prediction_train_val_c,opts.train.values);
        MABS(iter) = get_MABS(prediction_test_val_c,opts.test.values);
        MABS_train(iter) = get_MABS(prediction_train_val_c,opts.train.values);
        
        
    % check stopping
    if ((reschg < reschg_tol && ...
            itr_rank > minitr_reduce_rank) || relres < tol); break; end

    %if ((est_rank ==1||est_rank==0) && ((reschg < reschg_tol && ...
    %        itr_rank > minitr_reduce_rank) || relres < tol) ) ...
    %   || ((est_rank == 2)&&( relres < tol || (k==rank_max && itr_rank >= maxit)))
    %        break; 
    %end
    end
end %iter


if iprint == 1; fprintf('\n'); end
if data_tran; tX = X; X = Y'; Y = tX'; end
Out.obj = objv(1:iter);
Out.RR = RR(1:iter);
Out.iter = iter;
Out.rank = rk;
Out.relres = relres;
Out.reschg = reschg;
Out.time = time(1:iter);
if isfield(opts,'saveiterates') && opts.saveiterates == 1
    Out.Xhist=Out.Xhist(1:iter);
    Out.Yhist=Out.Yhist(1:iter);
else
    if m>n
        Out.Xhist=X;
        Out.Yhist=Y;
    else
        Out.Xhist=X;%Y';
        Out.Yhist=Y;%X';
    end
end
if recsys
   Out.RMSE      = RMSE(1:iter);
   Out.RMSE_train= RMSE_train(1:iter);
   Out.MABS      = MABS(1:iter);
   Out.MABS_train= MABS_train(1:iter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function init_inc()
        dr = L/(m+n-k)/k;  increment = .5*log10(m*n) - log10(rk);
        if min(m,n)/k > 1000; increment = increment + .25*exp(dr-1)/dr; end
    end

    function rank_estimator_adaptive()
        if est_rank == 1
            dR = abs(diag(R));       drops = dR(1:end-1)./dR(2:end);
            [dmx,imx] = max(drops);  rel_drp = (k-1)*dmx/(sum(drops)-dmx);
            %imx
            %bar(drops)
            if (rel_drp > rk_jump && itr_rank > minitr_reduce_rank) ...
                   || itr_rank > maxitr_reduce_rank; %bar(drops); pause;
                rk = max([imx, floor(0.1*k), rank_min]); 
                X = X(:,1:rk); Y = Y(1:rk,:);
                if Zfull %Z is full
                    Z = X*Y;  Res = data - Z(Known);   
                    Z(Known) = data + alf*Res; 
                else %Z = S + XY
                    Res = data - partXY(X',Y,Ik,Jk,L);
                    updateSval(S, (alf+1)*Res, L); 
                end
                res = norm(Res);  est_rank = 0; itr_rank = 0; 
                if iprint >= 2 %&& ch_rank>=1
                    fprintf('it: %5i, rel_drp: %3.2e, est_rank: %d,  Rank estimate changed from %i to %i\n',...
                        iter, rel_drp, est_rank, k,rk);
                end
            end
        elseif est_rank == 2 && reschg < 10*tol && rk < rank_max && itr_rank > 1 
        %elseif est_rank == 2 && reschg < tol/10  && itr_rank > 1 
            if rk < 50; rinc = rk_inc; else  rinc = 2*rk_inc; end
            rk = min(rk + rinc, rank_max); rkr_id = true;
            if rk > k;
                if save_res == 1
                    save(strcat('LM-Med-r', num2str(rk),'max', num2str(rank_max),  '.mat'), 'X', 'Y');
                end
                X = [X, zeros(m, rk-k)]; Y = [Y; zeros(rk-k,n)]; itr_rank = 0;
                if iprint >= 2 %&& ch_rank>=1
                    fprintf('it: %5i, reschg: %3.2e, Rank estimate changed from %i to %i\n',...
                        iter,reschg, k,rk);
                end
            end
        end
    end %rank

end %main
