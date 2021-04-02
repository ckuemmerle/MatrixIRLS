function vX = E2D_rankmanifold(gam,QU,QV)
% Implements adaption of Algorithm 6 (E2D_X^{M_p}) of 
% W. Huang, P.-A. Absil, K.A. Gullivan, "Intrinsic representation of 
% tangent vectors and vector transports on matrix manifolds", Numerische 
% Mathematik 136.2 (2017): 523-543.
% =========================================================================
% Based on code by Wen Huang, adapted by Christian Kuemmerle, Nov. 2020.
% =========================================================================
[d1,r] = size(QU.HHR);
d2 = size(QV.HHR,1);
[~,Gam3,Gam2] = get_Tk_matrices(gam,d1,d2,r);
K = HHRTeta(QU,Gam3);
K = K(r+1:end,:);
W = HHRTeta(QV,Gam2');
W = W(r+1:end,:)';

vX = zeros(r*(d1+d2-r),1);
vX(1:r^2) = gam(1:r^2);
vX((r^2+1):r*d1) = K(:);
vX(r*d1+1:end) = W(:);
end

