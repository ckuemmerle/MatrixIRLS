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
    gam =zeros(R*(d1+d2+R),1);
    setSval(sps_plc,y,m);
    UPhst_Yval = U'*sps_plc; 
    Phst_Yval = sps_plc*V;
    if strcmp(mode,'rangespace_smallsys')
    %%% Standard version for representations without U_{T_c} etc.:
    % y_T(1:R^2)=reshape(-UPhst_Yval*V,[R^2,1]);
    %%% Version for representation with U_{T_c} etc.:
        gam(1:R^2)=reshape(UPhst_Yval*V,[R^2,1]);
        gam((R^2+1):(R*(R+d2)))=reshape(UPhst_Yval,[R*d2,1]);
        gam((R*(R+d2)+1):end)=reshape(Phst_Yval,[d1*R,1]);
    elseif strcmp(mode,'tangspace')
        M1= UPhst_Yval*V;
        if increase_antisymmetricweights
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
            gam(1:R^2)              = reshape(M1,[R^2,1]);
            gam((R^2+1):(R*(R+d1))) = reshape(Phst_Yval-U*M1,[R*d1,1]);
            gam((R*(R+d1)+1):end)   = reshape(UPhst_Yval-M1*V',[R*d2,1]);
        end
    end
otherwise
    error('proj_tangspace_from_Oprange.m not yet implemented for this problem.')
end

end

