function [error_tail_rel,error_tail] =get_tail_frob_errors(Xr,X0,Phi,alg_name,k,frob_mode,verbose)
%get_tail_frob_errors This function computes the (relative) tail Frobenius
%discrepancies d_k(M) (which are not norms) of order k
% d_k(M) = \sqrt( \sum_{i=k+1}^d \sigma_i^2(M)) 
%        = \sqrt( \|M\|_F^2 - \sum_{i=1}^k \sigma_i^2(M))
%of the matrices M=Xr{l}{i} - X0.
if not(strcmp(frob_mode,'full'))
   disp('Discrepancies for partial matrices not implemented yet.')
end
thr_numinaccurate = 1e-4;


nr_algos=length(Xr);
error_tail       = cell(1,nr_algos);
error_tail_rel   = cell(1,nr_algos);
[d1,d2]=size(Phi);
meas = struct;
meas.Phi = Phi;
meas.Omega = find(Phi);
meas.Omega_c = setdiff([1:d1*d2],meas.Omega);
[meas.rowind,meas.colind]=find(Phi);
[meas.rowind_c,meas.colind_c] = ind2sub([d1,d2],meas.Omega_c);
meas.Phi_c = sparse(meas.rowind_c,meas.colind_c,ones(length(meas.Omega_c),1),d1,d2);

for l=1:nr_algos
    error_tail_rel{l}= zeros(1,length(Xr{l}));
    error_tail{l}    = zeros(1,length(Xr{l}));
    for i = 1:length(Xr{l})
        [error_tail_rel{l}(i),error_tail{l}(i)] = ...
        calc_tail_frob_error(Xr{l}{i},X0,meas,k,frob_mode,thr_numinaccurate);
    end
    if verbose
        disp(alg_name{l});
        disp(['Relative tail Frobenius error w.r.t. X0 (',frob_mode,'):  ',...
            num2str(error_tail_rel{l}(length(Xr{l})))]);
    end
    
    
end
end

function [error_tail_rel,error_fro_tail] = ...
    calc_tail_frob_error(Xr_c,X0,meas,k,frob_mode,thr_numinaccurate)

if isstruct(Xr_c)
    setSval(meas.Phi,Xr_c.res_range,length(Xr_c.res_range));
    if iscell(X0)
        U0=X0{1};
        V0=X0{2};
        if strcmp(frob_mode,'full')
            UU=U0'*Xr_c.U;
            VV=Xr_c.V'*V0;
            %%% Calculate the squared Frobenius distance:
            %%% ||X0-M||_F^2=||X0||_F^2+||M||_F^2 - 2<X0,M>
            normX0sq = trace((U0'*U0)*(V0'*V0));
            error_fro = norm_frob_compact(Xr_c,meas.Phi).^2 + normX0sq;
            error_fro = - 2* real(trace(VV*UU*Xr_c.Gam1))              + error_fro;
            error_fro2 = - 2* real(trace(UU*(Xr_c.Gam2*V0)));
            error_fro2 = - 2* real(trace(VV'*(Xr_c.Gam3'*U0)))          + error_fro2;
            error_fro2 = - 2* real(trace((V0'*meas.Phi')*U0))       + error_fro2;
            error_fro = error_fro + error_fro2;
            if error_fro/normX0sq < thr_numinaccurate || not(isreal(error_fro))
               error_fro= norm(get_densemat_from_compact(Xr_c,meas.Phi)-U0*V0','fro')^2;
            end
            M_handle = get_matrix_handle(Xr_c,meas.Phi,X0);
            [~, singval_M, ~]=bksvd_mod(M_handle, k, 200);
            error_fro_bulk = sum(diag(singval_M).^2);
            error_fro_tail = error_fro - error_fro_bulk;
            if error_fro_tail < 0
                error_fro_tail = 1e-16;
            end
        elseif strcmp(frob_mode,'Phi_comp') || strcmp(frob_mode,'Phi')
            if strcmp(frob_mode,'Phi')
%                     training_flag = 1;
%                     Xr_c_Phi = get_prediction_RecSys_IRLS({M},U,V, meas.rowind, ...
%                         meas.colind,[],[],training_flag);
%                     X0_Phi = partXY(U0', V0', meas.rowind, meas.colind, ...
%                     length(meas.colind))';
%                     error_fro = sumsqr(Xr_c_Phi - X0_Phi);
            elseif strcmp(frob_mode,'Phi_comp')
%                     training_flag = 0;
%                     Xr_c_Phi_c = get_prediction_RecSys_IRLS({M},U,V, meas.rowind_c, ...
%                         meas.colind_c,[],[],training_flag);
%                     X0_Phi_c = partXY(U0', V0', meas.rowind_c, meas.colind_c, ...
%                         length(meas.colind_c))';
%                     error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
            end
            error_fro_tail = 0;
        end
    else
        if strcmp(frob_mode,'full')
            %%% Calculate the squared Frobenius distance:
            %%% ||X0-M||_F^2=||X0||_F^2+||M||_F^2 - 2<X0,M>
            error_fro = norm_frob_compact(Xr_c,meas.Phi).^2;
            error_fro = trace(X0'*X0)                          + error_fro;
            error_fro = - 2* real(trace((Xr_c.V'*(X0'*Xr_c.U))*Xr_c.Gam1))       + error_fro;
            error_fro = - 2* real(trace((Xr_c.Gam2*X0')*Xr_c.U))            + error_fro;
            error_fro = - 2* real(trace((Xr_c.Gam3'*X0)*Xr_c.V))            + error_fro;
            error_fro = - 2* real(trace(meas.Phi'*X0))             + error_fro;
            if error_fro < thr_numinaccurate || not(isreal(error_fro))
               error_fro= norm(get_full_mat(Xr_c,meas.Phi)-X0,'fro')^2;
            end
            M_handle = get_matrix_handle(Xr_c,meas.Phi,X0);
            [~, singval_M, ~]=bksvd_mod(M_handle, k, 50);
            error_fro_bulk = sum(diag(singval_M).^2);
            error_fro_tail = error_fro - error_fro_bulk;
        elseif strcmp(frob_mode,'Phi')
%                 Xr_c_Phi = get_prediction_RecSys_IRLS({M},U,V, meas.rowind, ...
%                     meas.colind,[],[],training_flag);
%                 X0_Phi = X0(meas.Omega);
%                 error_fro = sumsqr(Xr_c_Phi - X0_Phi);
            error_fro_tail = 0;
        elseif strcmp(frob_mode,'Phi_comp')
%                 training_flag = 0;
%                 Xr_c_Phi_c = get_prediction_RecSys_IRLS({M},U,V, meas.rowind_c, ...
%                     meas.colind_c,[],[],training_flag);
%                 X0_Phi_c = X0(meas.Omega_c);
%                 error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
            error_fro_tail = 0;
        end
    end
%         if error_fro < 0
%             error_fro= norm(get_full_mat(Xr_c{1},Xr_c{2},Xr_c{3},meas.Phi)-U0*V0','fro'); 
%         else
    error_fro_tail =   sqrt(error_fro_tail);
%         end
elseif iscell(Xr_c)
    if size(Xr_c,2) == 2
        U=Xr_c{1};
        V=Xr_c{2};
        rX = size(U,2);
        if sum(sum(isnan(U)))
            U(isnan(U))=0;
        end
        if sum(sum(isnan(V)))
            V(isnan(V))=0;
        end
        
        if iscell(X0)
            U0=X0{1};
            r = size(U0,2); 
            V0=X0{2};
            if strcmp(frob_mode,'full')
%                 normX0sq = trace((U0'*U0)*(V0'*V0));
%                 error_fro = trace((U'*U)*(V'*V)) + normX0sq;
%                 error_fro = - 2* real(trace((U0'*U)*(V'*V0)))         + error_fro;
%                 error_fro = sqfrobnormfactors(U, V.', X0{1}, X0{2}.');
%                 if error_fro/normX0sq < thr_numinaccurate || not(isreal(error_fro))
%                     error_fro= norm(U*V'-X0{1}*X0{2}','fro')^2;
%                 end
                M_handle = get_matrix_handle({U,V},[],X0);
                [~, singval_M, ~]=bksvd_mod(M_handle, r+rX, 50);
                singval_M = diag(singval_M);
                last_singval_M=singval_M(k+1:(r+rX));
                %error_fro_bulk = sum(diag(singval_M).^2);
                %error_fro_tail = error_fro - error_fro_bulk;
                error_fro_tail = sum(last_singval_M.^2);
            elseif strcmp(frob_mode,'Phi')
%                 Xr_c_Phi = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
%                     meas.rowind,meas.colind,[],[]);
%                 X0_Phi = partXY(U0', V0', meas.rowind, meas.colind, ...
%                     length(meas.colind))';
%                 error_fro = sumsqr(Xr_c_Phi - X0_Phi);
                error_fro_tail = 0;
            elseif strcmp(frob_mode,'Phi_comp')
%                 Xr_c_Phi_c = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
%                     meas.rowind_c,meas.colind_c,[],[]);
%                 X0_Phi_c = partXY(U0', V0', meas.rowind_c, meas.colind_c, ...
%                     length(meas.colind_c))';
%                 error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
                error_fro_tail = 0;
            end
        else
            if strcmp(frob_mode,'full')
                error_fro = trace((U'*U)*(V'*V));
                error_fro = - 2* real(trace(V'*X0*U))                 + error_fro;
                error_fro = + norm(X0,'fro').^2                 + error_fro;
                M_handle = get_matrix_handle({U,V},[],X0);
                [~, singval_M, ~]=bksvd_mod(M_handle, k, 20);
                error_fro_bulk = sum(diag(singval_M).^2);
                error_fro_tail = error_fro - error_fro_bulk;
            elseif strcmp(frob_mode,'Phi')
%                 Xr_c_Phi = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
%                     meas.rowind,meas.colind,[],[]);
%                 X0_Phi = X0(meas.Omega);
%                 error_fro = sumsqr(Xr_c_Phi - X0_Phi);
                error_fro_tail = 0;
            elseif strcmp(frob_mode,'Phi_comp')
%                 Xr_c_Phi_c = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
%                     meas.rowind_c,meas.colind_c,[],[]);
%                 X0_Phi_c = X0(meas.Omega_c);
%                 error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
                error_fro_tail = 0;
            end
        end
        error_fro_tail =   sqrt(error_fro_tail);
    else
        error('Error in the format of the recovered matrices.')
        
    end
else
    if iscell(X0)
        U0=X0{1};
        V0=X0{2};
        X0_c = U0*V0';
    end
    if strcmp(frob_mode,'full')
        [d1,d2]=size(X0_c);
        delta_Xr_X0 = Xr_c-X0_c;
        error_fro = norm(delta_Xr_X0,'fro')^2;
        M_hdl  = @(x) (delta_Xr_X0)*x;
        Mt_hdl = @(y) (delta_Xr_X0)'*y;
        M_handle={M_hdl,Mt_hdl,d1,d2};
        [~, singval_M, ~]=bksvd_mod(M_handle, k, 20);
        error_fro_bulk = sum(diag(singval_M).^2);
        error_fro_tail = error_fro - error_fro_bulk;
        error_fro_tail =   sqrt(error_fro_tail);
    elseif strcmp(frob_mode,'Phi')
%         delta_Xr_X0 = Xr_c-X0;
%         error_fro = norm(delta_Xr_X0(meas.Omega),2);
        error_fro_tail = 0;
    elseif strcmp(frob_mode,'Phi_comp')
%         delta_Xr_X0 = Xr_c-X0;
%         error_fro = norm(delta_Xr_X0(meas.Omega_c),2);
        error_fro_tail = 0;
    end
end
if iscell(X0)
    if strcmp(frob_mode,'full')
        [d1,r]=size(X0{1});
        d2=size(X0{2},1);
        [UL,SL,VL]=svd(X0{1});
        [UR,SR,VR]=svd(X0{2});
        singval_X0=svds(SL(1:r,:)*VL'*VR*SR(1:r,:)',r);
%         X0_hdl  = @(x) X0{1}*(X0{2}'*x);
%         X0t_hdl = @(y) X0{2}*(X0{1}'*y);
%         X0_handle={X0_hdl,X0t_hdl,d1,d2};
%         [~, singval_X0, ~]=bksvd_mod(X0_handle, r, 50);
%         first_singval_X0 = singval_X0(1:k);
%         normX0sq = trace((U0'*U0)*(V0'*V0));
%         bulk_X0 = diag(singval_X0).^2;
        error_tail_rel = error_fro_tail./sqrt(sum(singval_X0(k+1:r).^2));
    elseif strcmp(frob_mode,'Phi')
%         error_fro_rel = error_fro./norm(get_prediction_RecSys_MatFac(U0,V0,eye(size(U0,2)),...
%                     meas.rowind,meas.colind,[],[]),2);
        error_tail_rel = 0;
    elseif strcmp(frob_mode,'Phi_comp')
%         error_fro_rel = error_fro./norm(get_prediction_RecSys_MatFac(U0,V0,eye(size(U0,2)),...
%                     meas.rowind_c,meas.colind_c,[],[]),2);
        error_tail_rel = 0;
    end
else
    if strcmp(frob_mode,'full')
        [d1,d2]=size(X0);
        X0_hdl  = @(x) X0*x;
        X0t_hdl = @(y) X0'*y;
        X0_handle={X0_hdl,X0t_hdl,d1,d2};
        [~, singval_X0, ~]=bksvd_mod(X0_handle, k, 50);
        singval_X0=diag(singval_X0);
        first_singval_X0=singval_X0(1:k);
        normX0sq = sum(sum(X0.^2));
%         bulk_X0 = diag(singval_X0).^2;
        error_tail_rel = error_fro_tail./(normX0sq - sum(first_singval_X0.^2));%(norm(X0,'fro')^2-sum(bulk_X0));
    elseif strcmp(frob_mode,'Phi')
%         error_fro_rel = error_fro./norm(X0(meas.Omega),2);
        error_tail_rel = 0;
    elseif strcmp(frob_mode,'Phi_comp')
%         error_fro_rel = error_fro./norm(X0(meas.Omega_c),2);
        error_tail_rel = 0;
    end
end

end

function M_handle = get_matrix_handle(M,M4_spsmat,X0)
if isstruct(M)
    d1 = size(M.U,1);
    d2 = size(M.V,1);
elseif iscell(M)
    U = M{1};
    V = M{2};
    d1 = size(U,1);
    d2 = size(V,1);
end
if iscell(M)
    if iscell(X0)
        M_hdl  = @(x) U*(V'*x) -X0{1}*(X0{2}'*x);
        Mt_hdl = @(y) V*(U'*y)-X0{2}*(X0{1}'*y);
    else
        M_hdl  = @(x) U*(V'*x) -X0*x;
        Mt_hdl = @(y) V*(U'*y)-X0'*y;
    end
else
    if iscell(X0)
        M_hdl  = @(x) (M.U*(M.Gam1*(M.V'*x) +(M.Gam2*x)) + M.Gam3*(M.V'*x) +M4_spsmat*x...
            -X0{1}*(X0{2}'*x));
        Mt_hdl = @(y) (M.V*(M.Gam1'*(M.U'*y)+M.Gam3'*y) + M.Gam2'*(M.U'*y) +M4_spsmat'*y...
            -X0{2}*(X0{1}'*y));
    else
        M_hdl  = @(x) (M.U*(M.Gam1*(M.V'*x) +(M.Gam2*x)) + M.Gam3*(M.V'*x) +M4_spsmat*x...
            -X0*x);
        Mt_hdl = @(y) (M.V*(M.Gam1'*(M.U'*y)+M.Gam3'*y) + M.Gam2'*(M.U'*y) +M4_spsmat'*y...
            -X0'*y);
    end
end
M_handle={M_hdl,Mt_hdl,d1,d2};
end