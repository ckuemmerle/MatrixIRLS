function [X_Tk_1,X_Tk_2,X_Tk_3] = get_Tk_matrices(gam,d1,d2,r)
%get_Tk_matrices Summary of this function goes here
%   Detailed explanation goes here
if not(length(gam) == r*(d1+d2+r))
    error('Error in the dimensionality.')
end
X_Tk_1=reshape(gam(1:r^2),[r,r]);
X_Tk_2=reshape(gam((r^2+1):(r*(r+d1))),[d1,r]);
X_Tk_3=reshape(gam((r*(r+d1)+1):end),[r,d2]);
end

