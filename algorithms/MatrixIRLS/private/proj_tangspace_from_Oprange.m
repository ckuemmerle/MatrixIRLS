function gam = proj_tangspace_from_Oprange(y,problem,varargin)
%proj_tangentspace_from_Omega This function projects a m-sparse matrix y
%with values in Omega
%into the tangent space
% T:= {U*Z1' + Z2*V' ; Z1 \in \R^(d1xr), Z2 \in \R^(d2xr)}. The map can be
% written as
%P_T(P_Omega'(y)) = U*U'*P_Omega'(y) + *P_Omega'(y)*V*V' - U*U'**P_Omega'(y)*V*V'.
switch problem
case 'MatrixCompletion'
    U       = varargin{1};
    V       = varargin{2};
    sps_plc = varargin{3};
    mode    = varargin{4};
    if nargin == 7
       increase_antisymmetricweights = varargin{5}; 
    else
       increase_antisymmetricweights = 0;
    end
    
    m=length(y);
    [d1,R]=size(U);
    d2    =size(V,1);
    D=max(d1,d2);
    setSval(sps_plc,y,m);
%     vals=nonzeros(sps_plc);
%     UPhst_Yval_tst = multfullsparse(U', vals, sps_plc);
%     Phst_Yval_tst = multsparsefull(vals, V, sps_plc); (this is not
%     faster)
    if strcmp(mode,'rangespace_smallsys')
        UPhst_Yval = U'*sps_plc;
        Phst_Yval = sps_plc*V;
    %%% Standard version for representations without U_{T_c} etc.:
    % y_T(1:R^2)=reshape(-UPhst_Yval*V,[R^2,1]);
    %%% Version for representation with U_{T_c} etc.:
        gam =zeros(R*(d1+d2+R),1);
        gam(1:R^2)=reshape(UPhst_Yval*V,[R^2,1]);
        gam((R^2+1):(R*(R+d2)))=reshape(UPhst_Yval,[R*d2,1]);
        gam((R*(R+d2)+1):end)=reshape(Phst_Yval,[d1*R,1]);
    elseif strcmp(mode,'tangspace')
        if increase_antisymmetricweights
            UPhst_Yval = U'*sps_plc;
            Phst_Yval = sps_plc*V;
            M1= UPhst_Yval*V;
            gam =zeros(R*(d1+d2+R),1);
            M1S = (M1+M1')./2;
            M1T = (M1-M1')./2;
            M1_upper  = triu(M1S);
            M1_lower = tril(M1T,-1);
            M1  = M1_upper+M1_lower;
            if d1 == D
                M2=(Phst_Yval-U*M1+[UPhst_Yval'-V*M1' ;zeros(d1-d2,d1)])./2;
                tmp=Phst_Yval'-M1'*U';
                M3=(UPhst_Yval-M1*V'-tmp(:,1:d2))./2;
                gam(1:R^2)              = reshape(M1,[R^2,1]);
                gam((R^2+1):(R*(R+d1))) = reshape(M2,[R*d1,1]);
                gam((R*(R+d1)+1):end)   = reshape(M3,[R*d2,1]);
            else
               error('To be implemented.') 
            end
        else
            gam = matrixspace_to_tangspace(sps_plc,U,V);
%             gam(1:R^2)              = reshape(M1,[R^2,1]);
%             gam((R^2+1):(R*(R+d1))) = reshape(Phst_Yval-U*M1,[R*d1,1]);
%             gam((R*(R+d1)+1):end)   = reshape(UPhst_Yval-M1*V',[R*d2,1]);
        end
    end
    
case 'RobustPCA'
    U       = varargin{1};
    V       = varargin{2};
%     [d1,R]=size(U);
%     d2    =size(V,1);
    [d1,d2] = size(y);
    if isequal(size(y),[d1,d2])
        Y=y;
    else
        Y=reshape(y,[d1,d2]);
    end
    gam = matrixspace_to_tangspace(Y,U,V);
%     R = size(U,2);
%     gam =zeros(R*(d1+d2+R),1);
%     UPhst_Yval= U'*Y;
%     Phst_Yval = Y*V;
%     M1= UPhst_Yval*V;
%     gam(1:R^2)              = reshape(M1,[R^2,1]);
%     gam((R^2+1):(R*(R+d1))) = reshape(Phst_Yval-U*M1,[R*d1,1]);
%     gam((R*(R+d1)+1):end)   = reshape(UPhst_Yval-M1*V',[R*d2,1]);
case 'PhaseRetrieval'
    U  = varargin{1};
    At = varargin{2};
    AU = varargin{3};
    [~,R] = size(AU);
	n = size(U,1);
    
    yAU=y.*AU;
    X1=AU'*(yAU);
    handle_At=functions(At);
    if contains(handle_At.function,'transposeOperator')
        At_op= @(X) cell2mat(cellfun(At,num2cell(X,1),'UniformOutput',false));
        AtyAU=conj(At_op(conj(yAU)));
        X2=sqrt(2).*(AtyAU-U*(AU'*yAU));
    else
        X2=sqrt(2).*(conj(At(conj(yAU)))-U*(AU'*yAU));
    end
%     X2=sqrt(2).*(conj(At(conj(yAU)))-U*(AU'*yAU)); %At_op(conj(yAU)
    gam = zeros(R*(n+R),1);
    gam(1:R^2)=reshape(X1,[R^2,1]);
    gam((R^2+1):(R*(R+n)))=reshape(X2,[R*n,1]);
    
otherwise
    error('proj_tangspace_from_Oprange.m not yet implemented for this problem.')
end

end

