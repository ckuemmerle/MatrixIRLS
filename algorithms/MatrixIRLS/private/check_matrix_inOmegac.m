function norm1 = check_matrix_inOmegac(M,U,V,rowind,colind,sparse_placemat)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
Mnew=proj_Omega(M,U,V,rowind,colind);
X = get_full_mat(Mnew,U,V,sparse_placemat);
disp('Energy in the Omega subspace:')
norm1=norm(X,'fro');

end

