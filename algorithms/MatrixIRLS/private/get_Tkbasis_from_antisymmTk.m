function gam_Tk = get_Tkbasis_from_antisymmTk(gam,d1,d2,r)
%get_Tkbasis_from_antisymmTk Summary of this function goes here
%   Detailed explanation goes here
gam_Tk=gam;
[M1,M2,M3] = get_Tk_matrices(gam,d1,d2,r);

M1S_u = triu(M1);
M1S_l = triu(M1,1)';
M1S   = M1S_u+M1S_l;
M1T_l = tril(M1,-1);
M1T_u = -M1T_l';    
M1T   = M1T_u+M1T_l;
Z1    = M1S+M1T;

if d1 == max(d1,d2)
    Z2    = M2+[M3';zeros(d1-d2,r)];
    tmp=M2;
    Z3    = -M3+tmp(1:d2,:)';
else
    error('To be implemented.')
end

gam_Tk(1:r^2)              = reshape(Z1,[r^2,1]);
gam_Tk((r^2+1):(r*(r+d1))) = reshape(Z2,[r*d1,1]);
gam_Tk((r*(r+d1)+1):end)   = reshape(Z3,[r*d2,1]);

end

