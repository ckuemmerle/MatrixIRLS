function gam = matrixspace_to_tangspace(X,U,V)
%matrixspace_to_tangspace Corresponds to operator "P_{T_k}^*".
[d1,r]=size(U);
d2    =size(V,1);
gam =zeros(r*(d1+d2+r),1);
UtX = U'*X; 
XV = X*V;
M1= UtX*V;
gam(1:r^2)              = reshape(M1,[r^2,1]);
gam((r^2+1):(r*(r+d1))) = reshape((XV-U*M1),[r*d1,1]);
gam((r*(r+d1)+1):end)   = reshape((UtX-M1*V'),[r*d2,1]);
end