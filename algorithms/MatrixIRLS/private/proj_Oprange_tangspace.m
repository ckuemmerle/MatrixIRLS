function y = proj_Oprange_tangspace(gam,problem,varargin)
%proj_Omega_MC_smallsys This function projects a rank-2R (spaces (U,V))
%matrix onto its entries with indices in Omega
%       Input:      gamma = (R*(d1+d2-R))x1 vector. (for matrix completion)
if strcmp(problem,'MatrixCompletion')
    U      = varargin{1};
    V      = varargin{2};
    rowind = varargin{3};
    colind = varargin{4};
    mode   = varargin{5};
    if nargin == 8
       increase_antisymmetricweights = varargin{6}; 
    else
       increase_antisymmetricweights = 0;
    end
    
    [d1,R]  =size(U);
    d2      =size(V,1);
    D=max(d1,d2);
    m=length(rowind);

    %%% Standard version for representations without U_{T_c} etc.:
    % Z2 = reshape(gam((R^2+1):(R*(d2+R))),[R,d2]);
    % Z3 = reshape(gam((R*(d2+R)+1):(R*(d2+d1+R))),[d1,R]);
    % if isreal(gam)
    %     gam_Omega = partXY((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY(U.',Z2,rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY((Z3).',V',rowind,colind,m)';
    % else
    %     gam_Omega = partXY_cmplx((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY_cmplx(U.',Z2,rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY_cmplx((Z3).',V',rowind,colind,m)';
    % end
    %%% Version where a (R x R) block is set 0 -> problem later on:
    % Z2 = [zeros(R,R),reshape(gam((R^2+1):(R*(d2))),[R,d2-R])];
    % Z3 = [zeros(R,R);reshape(gam((R*d2+1):(R*(d2+d1-R))),[d1-R,R])];
    %%% Version without 0 setting, but still representation with U_{T_c}:
    switch mode
        case 'rangespace_smallsys'
            M2 = reshape(gam((R^2+1):(R*(d2+R))),[R,d2]);
            M3 = reshape(gam((R*(d2+R)+1):(R*(d2+d1+R))),[d1,R]);
            Z2V = M2*V;
            UZ3 = U'*M3;

            if isreal(gam)
                y = partXY((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
                y = y+partXY(U.',M2-Z2V*V',rowind,colind,m)';
                %   gam_Omega+partXY(U.',reshape(gam((R^2+1):(R*(d2))),[R,d2]),rowind,colind,m)';
                y = y+partXY((M3-U*UZ3).',V',rowind,colind,m)';
                %   gam_Omega = gam_Omega+partXY((reshape(gam((R*(R+d2)+1):end),[d1,R])).',V',rowind,colind,m)';
            else
                y = partXY_cmplx((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
                y = y+partXY_cmplx(U.',M2-Z2V*V',rowind,colind,m)';
                y = y+partXY_cmplx((M3-U*UZ3).',V',rowind,colind,m)';
            end
        case 'tangspace'
            [M1,M2,M3] = get_Tk_matrices(gam,d1,d2,R);
%             M1=reshape(gam(1:R^2),[R,R]);
%             M2=reshape(gam((R^2+1):(R*(d1+R))),[d1,R]);
%             M3=reshape(gam((R*(d1+R)+1):(R*(d2+d1+R))),[R,d2]);
            if increase_antisymmetricweights
                M1S_u = triu(M1);
                M1S_l = triu(M1,1)';
                M1S   = M1S_u+M1S_l;
                M1T_l = triu(M1,-1);
                M1T_u = -M1T_l';
                M1T   = M1T_u+M1T_l;
                Z1    = M1S+M1T;
                
                if d1 == D
                    Z2    = M2+[M3',zeros(d1,d1-d2)];
                    tmp=M2;
                    Z3    = M3+tmp(1:d2,:)';
                else
                    error('To be implemented.')
                end
                M1 = Z1;
                M2 = Z2;
                M3 = Z3;
            end
            if isreal(gam)
%                 rowindit = uint32(rowind);
%                 colindit = uint32(colind);
%                 tmp = spmaskmult((U*M1+M2), V', rowindit, colindit); this
%                 is not faster...
                y = partXY((U*M1+M2)',V',rowind,colind,m)';
                y = y+partXY(U.',M3,rowind,colind,m)';
    %             y = partXY((U*reshape(gam(1:R^2),[R,R])+reshape(gam((R^2+1):(R*(d1+R))),[d1,R])).',V',rowind,colind,m)';
    %             y = y+partXY(U.',reshape(gam((R*(d1+R)+1):(R*(d2+d1+R))),[R,d2]),rowind,colind,m)';
            else
                y = partXY_cmplx((U*M1+M2).',V',rowind,colind,m)';
                y = y+partXY_cmplx(U.',M3,rowind,colind,m)';
            end
    end
    
elseif strcmp(problem,'PhaseRetrieval')
    U       = varargin{1};
    A       = varargin{2};
    AU      = varargin{3};
    [~,R] = size(AU);
    n   = size(U,1);
%   X1= reshape(gam(1:R^2),[R,R]);
%   X2= related to reshape(gam((R^2+1):end),[n,R]);
    handle_A=functions(A);
    if contains(handle_A.function,'fourierMeasurementOperator')
        Aop= @(X) cell2mat(cellfun(A,num2cell(X,1),'UniformOutput',false));
        tmp=Aop(reshape(conj(gam((R^2+1):end)),[n,R]));
    else
        tmp=A(reshape(conj(gam((R^2+1):end)),[n,R]));
    end
    y = sum((AU*reshape(gam(1:R^2),[R,R])).*conj(AU),2);
    y = y + sqrt(2).*real(sum(AU.*tmp,2));
%     y = y + sqrt(2).*real(sum(AU.*(Aop(reshape(conj(gam((R^2+1):end)),[n,R]))),2));
else
    error('proj_Oprange_tangspace.m not yet implemented for this problem.')
end
    
end

