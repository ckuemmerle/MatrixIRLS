function gam_T = proj_tangspace(gam,problem,U0,V0,U,V)
%proj_tangspace This function projects the matrices from the tangent space
%T0 into the tangent space T

Utild_U = U'*U0;
if strcmp(problem,'MatrixCompletion')
    [d1,r_T0]=size(U0);
    [~,r_T]=size(U);
    d2    =size(V0,1);
    V_Vtild = V0'*V;
    gam_T=zeros(r_T*(d1+d2+r_T),1); %gam;
    X1= reshape(gam(1:r_T0^2),[r_T0,r_T0]);
    X2= reshape(gam((r_T0^2+1):(r_T0*(d1+r_T0))),[d1,r_T0]);
    X3= reshape(gam((r_T0*(d1+r_T0)+1):(r_T0*(d2+d1+r_T0))),[r_T0,d2]);
    gam_T_1=Utild_U*(X1*V_Vtild+X3*V)+U'*X2*V_Vtild;
    gam_T_2=(U0-U*Utild_U)*(X3*V)+(U0*X1+X2-U*Utild_U*X1-U*(U'*X2))*V_Vtild;
    gam_T_3=Utild_U*(X1*V0'+X3-X1*V_Vtild*V'-(X3*V)*V')+...
        U'*X2*(V0'-V_Vtild*V');
    gam_T(1:r_T^2)                = reshape(gam_T_1,[r_T^2,1]);
    gam_T((r_T^2+1):(r_T*(r_T+d1)))   = reshape(gam_T_2,[r_T*d1,1]);
    gam_T((r_T*(r_T+d1)+1):end)     = reshape(gam_T_3,[r_T*d2,1]);
else
    error('proj_tangspace.m not yet implemented for this problem.')
end

end