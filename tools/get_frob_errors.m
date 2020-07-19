function [error_fro_rel,error_fro] =get_frob_errors(Xr,X0,Phi,alg_name,frob_mode,verbose)

thr_numinaccurate = 1e-3;

nr_algos=length(Xr);
error_fro       = cell(1,nr_algos);
error_fro_rel   = cell(1,nr_algos);
[d1,d2]=size(Phi);
meas = struct;
meas.Phi = Phi;
meas.Omega = find(Phi);
meas.Omega_c = setdiff([1:d1*d2],meas.Omega);
[meas.rowind,meas.colind]=find(Phi);
[meas.rowind_c,meas.colind_c] = ind2sub([d1,d2],meas.Omega_c);
meas.Phi_c = sparse(meas.rowind_c,meas.colind_c,ones(length(meas.Omega_c),1),d1,d2);

for l=1:nr_algos
    error_fro_rel{l}= zeros(1,length(Xr{l}));
    error_fro{l}    = zeros(1,length(Xr{l}));
    for i = 1:length(Xr{l})
        [error_fro_rel{l}(i),error_fro{l}(i)] = ...
        calc_frob_error(Xr{l}{i},X0,meas,frob_mode,thr_numinaccurate);
    end
    if verbose
        disp(alg_name{l});
        disp(['Relative Frob. error to X0 (',frob_mode,'):  ',num2str(error_fro_rel{l}(length(Xr{l})))]);
    end
    
    
end
end

function [error_fro_rel,error_fro] = ...
    calc_frob_error(Xr_c,X0,meas,frob_mode,thr_numinaccurate)

if isstruct(Xr_c)
        setSval(meas.Phi,Xr_c.res_range,length(Xr_c.res_range));
        if iscell(X0)
            U0=X0{1};
            V0=X0{2};
            if strcmp(frob_mode,'full')
                norm_X0 = sqrt(real(trace((U0'*U0)*(V0'*V0))));
                
                UU=U0'*Xr_c.U;
                VV=Xr_c.V'*V0;
                %%% Calculate the squared Frobenius distance:
                %%% ||X0-M||_F^2=||X0||_F^2+||M||_F^2 - 2<X0,M>
                error_fro = norm_frob_compact(Xr_c,meas.Phi).^2;
                error_fro = trace((U0'*U0)*(V0'*V0))            + error_fro;
                error_fro = - 2* real(trace(VV*UU*Xr_c.Gam1))              + error_fro;
                error_fro = - 2* real(trace(UU*(Xr_c.Gam2*V0)))            + error_fro;
                error_fro = - 2* real(trace(VV'*(Xr_c.Gam3'*U0)))          + error_fro;
                error_fro = - 2* real(trace((V0'*meas.Phi')*U0))       + error_fro;
                if error_fro/norm_X0 < thr_numinaccurate || not(isreal(error_fro))
%                     disp(['Recalculate diff_c for it = ',num2str(l)])
                   error_fro= norm(get_densemat_from_compact(Xr_c,meas.Phi)-U0*V0','fro')^2;
                end
            elseif strcmp(frob_mode,'Phi_comp') || strcmp(frob_mode,'Phi')
                if strcmp(frob_mode,'Phi')
                    training_flag = 1;
                    Xr_c_Phi = get_prediction_RecSys_IRLS(Xr_c, meas.rowind, ...
                        meas.colind,[],[],training_flag);
                    X0_Phi = partXY(U0', V0', meas.rowind, meas.colind, ...
                    length(meas.colind))';
                    error_fro = sumsqr(Xr_c_Phi - X0_Phi);
                elseif strcmp(frob_mode,'Phi_comp')
                    training_flag = 0;
                    Xr_c_Phi_c = get_prediction_RecSys_IRLS(Xr_c, meas.rowind_c, ...
                        meas.colind_c,[],[],training_flag);
                    X0_Phi_c = partXY(U0', V0', meas.rowind_c, meas.colind_c, ...
                        length(meas.colind_c))';
                    error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
                end
            end
        else
            if strcmp(frob_mode,'full')
                norm_X0 = norm(X0,'fro');
                %%% Calculate the squared Frobenius distance:
                %%% ||X0-M||_F^2=||X0||_F^2+||M||_F^2 - 2<X0,M>
                error_fro = norm_frob_compactformat(Xr_c,meas.Phi).^2;
                error_fro = trace(X0'*X0)                          + error_fro;
                error_fro = - 2* real(trace((V'*(X0'*U))*Xr_c.Gam1))       + error_fro;
                error_fro = - 2* real(trace((Xr_c.Gam2*X0')*U))            + error_fro;
                error_fro = - 2* real(trace((Xr_c.Gam3'*X0)*V))            + error_fro;
                error_fro = - 2* real(trace(meas.Phi'*X0))             + error_fro;
            elseif strcmp(frob_mode,'Phi')
                training_flag = 1;
                Xr_c_Phi = get_prediction_RecSys_IRLS(Xr_c, meas.rowind, ...
                    meas.colind,[],[],training_flag);
                X0_Phi = X0(meas.Omega);
                error_fro = sumsqr(Xr_c_Phi - X0_Phi);
            elseif strcmp(frob_mode,'Phi_comp')
                training_flag = 0;
                Xr_c_Phi_c = get_prediction_RecSys_IRLS(Xr_c, meas.rowind_c, ...
                    meas.colind_c,[],[],training_flag);
                X0_Phi_c = X0(meas.Omega_c);
                error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
            end
        end
        if error_fro < 0
            error_fro= norm(get_densemat_from_compact(Xr_c, meas.Phi)...
                       - U0*V0','fro'); 
        else
            error_fro =   sqrt(error_fro);
        end
elseif iscell(Xr_c)
    if size(Xr_c,2) == 2
        U=Xr_c{1};
        V=Xr_c{2};
        
        if iscell(X0)
            U0=X0{1};
            V0=X0{2};
            if strcmp(frob_mode,'full')
                norm_X0 = sqrt(real(trace((U0'*U0)*(V0'*V0))));
                error_fro = trace((U'*U)*(V'*V));
                error_fro = - 2* real(trace((U0'*U)*(V'*V0)))         + error_fro;
                error_fro = + trace((U0'*U0)*(V0'*V0))          + error_fro;
                if error_fro < thr_numinaccurate || not(isreal(error_fro))
                    error_fro= norm(U*V'-X0{1}*X0{2}','fro')^2;
                end
            elseif strcmp(frob_mode,'Phi')
                Xr_c_Phi = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
                    meas.rowind,meas.colind,[],[]);
                X0_Phi = partXY(U0', V0', meas.rowind, meas.colind, ...
                    length(meas.colind))';
                error_fro = sumsqr(Xr_c_Phi - X0_Phi);
            elseif strcmp(frob_mode,'Phi_comp')
                Xr_c_Phi_c = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
                    meas.rowind_c,meas.colind_c,[],[]);
                X0_Phi_c = partXY(U0', V0', meas.rowind_c, meas.colind_c, ...
                    length(meas.colind_c))';
                error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
            end
        else
            if strcmp(frob_mode,'full')
                norm_X0 = norm(X0,'fro');
                error_fro = trace((U'*U)*(V'*V));
                error_fro = - 2* real(trace(V'*X0*U))                 + error_fro;
                error_fro = + norm(X0,'fro').^2                 + error_fro;
            elseif strcmp(frob_mode,'Phi')
                Xr_c_Phi = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
                    meas.rowind,meas.colind,[],[]);
                X0_Phi = X0(meas.Omega);
                error_fro = sumsqr(Xr_c_Phi - X0_Phi);
            elseif strcmp(frob_mode,'Phi_comp')
                Xr_c_Phi_c = get_prediction_RecSys_MatFac(U,V,eye(size(U,2)),...
                    meas.rowind_c,meas.colind_c,[],[]);
                X0_Phi_c = X0(meas.Omega_c);
                error_fro = sumsqr(Xr_c_Phi_c - X0_Phi_c);
            end
        end
        error_fro =   sqrt(error_fro);
    else
        error('Error in the format of the recovered matrices.')
        
    end
else
    if iscell(X0)
        U0=X0{1};
        V0=X0{2};
        X0 = U0*V0';
    end
    if strcmp(frob_mode,'full') 
        norm_X0 = norm(X0,'fro');
        error_fro = norm(Xr_c-X0,'fro');
    elseif strcmp(frob_mode,'Phi')
        delta_Xr_X0 = Xr_c-X0;
        error_fro = norm(delta_Xr_X0(meas.Omega),2);
    elseif strcmp(frob_mode,'Phi_comp')
        delta_Xr_X0 = Xr_c-X0;
        error_fro = norm(delta_Xr_X0(meas.Omega_c),2);
    end
end
error_fro=real(error_fro);
if iscell(X0)
    if strcmp(frob_mode,'full')
        error_fro_rel = error_fro./norm_X0;
    elseif strcmp(frob_mode,'Phi')
        error_fro_rel = error_fro./norm(get_prediction_RecSys_MatFac(U0,V0,eye(size(U0,2)),...
                    meas.rowind,meas.colind,[],[]),2);
    elseif strcmp(frob_mode,'Phi_comp')
        error_fro_rel = error_fro./norm(get_prediction_RecSys_MatFac(U0,V0,eye(size(U0,2)),...
                    meas.rowind_c,meas.colind_c,[],[]),2);
    end
else
    if strcmp(frob_mode,'full')
        error_fro_rel = error_fro./norm_X0;
    elseif strcmp(frob_mode,'Phi')
        error_fro_rel = error_fro./norm(X0(meas.Omega),2);
    elseif strcmp(frob_mode,'Phi_comp')
        error_fro_rel = error_fro./norm(X0(meas.Omega_c),2);
    end
end

end