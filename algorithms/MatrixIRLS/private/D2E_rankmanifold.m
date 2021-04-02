function gam = D2E_rankmanifold(vX,QU,QV)
% Implements adaption of Algorithm 7 (D2E_X^{M_p}) of 
% W. Huang, P.-A. Absil, K.A. Gullivan, "Intrinsic representation of 
% tangent vectors and vector transports on matrix manifolds", Numerische 
% Mathematik 136.2 (2017): 523-543.
% =========================================================================
% Based on code by Wen Huang, adapted by Christian Kuemmerle, Nov. 2020.
% =========================================================================
[d1,r] = size(QU.HHR);
d2 = size(QV.HHR,1);
gam = zeros(r*(d1+d2+r),1);

% [~,Gam3,Gam2] = get_Tk_matrices(gam,d1,d2,r);

Klong = [zeros(r,r);reshape(vX((r^2+1):r*d1),[d1-r,r])];
Kback = HHReta(QU,Klong);

Wlong = [zeros(r,r);reshape(vX(r*d1+1:end),[r,d2-r])'];
Wback = HHReta(QV,Wlong);
% 
% X_Tk_2=reshape(gam((r^2+1):(r*(r+d1))),[d1,r]);
% X_Tk_3=reshape(gam((r*(r+d1)+1):end),[r,d2]);
gam(1:r^2) = vX(1:r^2);
gam(r^2+1:r*(r+d1)) = reshape(Kback,[r*d1,1]);
gam(r*(r+d1)+1:end) = reshape(Wback',[r*d2,1]);

% vX = zeros(r*(d1+d2-r),1);
% vX(1:r^2) = gam(1:r^2);
% vX((r^2+1):r*d1) = K(:);
% vX(r*d1+1:end) = W(:);
end

