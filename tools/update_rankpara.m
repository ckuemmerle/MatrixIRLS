function [r_c,sing,U_X,V_X]=update_rankpara(XXX,sing,U_X,V_X,eps,R,r_c,tracking)
%update_rankpara This functions determines how many singular values
%are larger than eps, and updates the rank r_c accordingly. It further
%outputs the first r_c singular values (sing) and first r_c singular vectors
%(U_X and V_X).
%
%       Input:  XXX = (1 x 4) cell with function handles of X and X', 
%                     respectively, and the dimensions d1 and d2.
%                 R = maximal number of singular values to be computed.
%               r_c = Number of singular values of last iterate larger than
%                     the eps of last iterate.
%                     
tol_precision = 1e-10;

d1=size(U_X,1);
d2=size(V_X,1);
d=min(d1,d2);
r_SVD=max(r_c,1);
r_correct = 0;
while not(r_correct)
    if length(find(sing>eps+tol_precision)) <= r_SVD
        r_c = length(find(sing>eps+tol_precision));
        r_correct = 1;
    else
        r_SVD = min(2*r_SVD,R);
        if ~tracking
            [U_X, singval_mat_c, V_X]=bksvd_mod(XXX, min(d,r_SVD+1), 20);
            sing = diag(singval_mat_c);
        end
        if length(find(sing>eps+tol_precision)) <= r_SVD
            r_c = length(find(sing>eps+tol_precision));
            r_correct = 1;
        end
%             else
%                 error('Implement loop here to calculate even more singular values.')
%             end
    end
    if r_SVD == R
        r_c=min(r_c,r_SVD);
        r_correct = 1;
    end
end
if r_c == 0
    sing = [];%eps;
    U_X = zeros(d1,r_c);
    V_X = zeros(d2,r_c);
else
    sing = sing(1:r_c);
    U_X = U_X(:,1:r_c);
    V_X = V_X(:,1:r_c);
end
end

