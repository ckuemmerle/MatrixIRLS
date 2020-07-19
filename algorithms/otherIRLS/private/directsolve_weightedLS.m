function X_matrix = directsolve_weightedLS(PhWiPh_star,Y,W_l_inv,W_r_inv,Omega,variant,varargin)
%it_min_l2_w_fast_MC solves the least squares problem of the IRLS, adapted
%to the matrix completion setting.

%%% PhWiPh_star : (m x m) matrix
%%% Y           : (m x 1) vector
%%% W_l_inv     : (d1 x d1) matrix
%%% W_r_inv     : (d2 x d2) matrix
%%% Omega       : (m x 1) vector containing the linear indices of the
%%%               revealed entries

%%% X_matrix    : (d1 x d2) matrix
if strcmp(variant,'arithmetic')
   U = W_l_inv;
   V = W_r_inv;
   H = varargin{1};
   tmp3=PhWiPh_star\Y;
   tmp4= zeros(size(U,1),size(V,1));
   tmp4(Omega(:))=tmp3(:);
else
%    tmp3=PhWiPh_star\sparse(Y);
   tmp3=PhWiPh_star\Y;
   tmp4 = sparse(size(W_l_inv,1)*size(W_r_inv,1),1);
   tmp4(Omega(:))= tmp3(:);
   tmp4_resh=reshape(tmp4,[size(W_l_inv,1),size(W_r_inv,1)]); 
end


if strcmp(variant,'harmonic')
    X_matrix = (0.5).*(W_l_inv*tmp4_resh+tmp4_resh*W_r_inv);
elseif strcmp(variant,'left')
    X_matrix = W_l_inv*tmp4_resh;
elseif strcmp(variant,'right')
    X_matrix = tmp4_resh*W_r_inv;
elseif strcmp(variant,'arithmetic')
    X_matrix = U*(H.*(U'*tmp4*V))*V';
end
end
