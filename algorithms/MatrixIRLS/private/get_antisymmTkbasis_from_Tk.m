function gam_Tkanti = get_antisymmTkbasis_from_Tk(gam_Tk,d1,d2,r)
%get_Tkbasis_from_antisymmTk Summary of this function goes here
%   Detailed explanation goes here
gam_Tkanti=gam_Tk;

[Z1,Z2,Z3] = get_Tk_matrices(gam_Tk,d1,d2,r);

Z1S = (Z1+Z1')./2;
Z1T = (Z1-Z1')./2;
M1_upper  = triu(Z1S);
M1_lower = tril(Z1T,-1);
M1  = M1_upper+M1_lower;

if d1 == max(d1,d2)
    M2=(Z2+[Z3' ;zeros(d1-d2,r)])./2;
    tmp=Z2';
    M3=(tmp(:,1:d2)-Z3)./2;
    gam_Tkanti(1:r^2)              = reshape(M1,[r^2,1]);
    gam_Tkanti((r^2+1):(r*(r+d1))) = reshape(M2,[r*d1,1]);
    gam_Tkanti((r*(r+d1)+1):end)   = reshape(M3,[r*d2,1]);
else
   error('To be implemented.') 
end

end

