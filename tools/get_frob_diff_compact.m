function err_fro = get_frob_diff_compact(M,N,sps_plc)%(M,U0,V0,N,U,V,sps_plc)
% Calculates the Frobenius distance ||M-N||_F 
% between two matrices M and N that are given as in the Appendix of 
% [1] C. Kuemmerle, C. M. Verdun, "Escaping Saddle Points in 
% Ill-Conditioned Matrix Completion with a Scalable Second Order Method", 
% ICML 2020 Workshop "Beyond First Order Methods in ML Systems".
% =========================================================================
% Parameters
% ----------
% M:    struct. Represents a matrix M_mat, is given such that,
%       M.Gam1 = \tilde{Gamma}_1  [(r_M x r_M) matrix]
%       M.Gam2 = \tilde{Gamma}_2  [(r_M x d_2) matrix]
%       M.Gam3 = \tilde{Gamma}_3  [(d_1 x r_M) matrix]
%       M.res_range = r_{k+1}     [(m x 1) vector].
%       M.U    = U^{(k)}          [(d_1 x r_M) matrix w/
%                                   orthonormal columns]
%       M.V    = V^{(k)}          [(d_2 x r_M) matrix w/
%                                   orthonormal columns],
%       such that its matrix form can be written as
%       M_mat =  U^{(k)}*(\tilde{Gamma}_1*V^{(k)'}+\tilde{Gamma}_2)
%           + \tilde{Gamma}_3*V^{(k)'} + P_{\Omega}^*(r_{k+1}).
% N:    struct. See parameter M.
% sps_plc: (d_1 x d_2) sparse matrix with support of size m. Non-zeros 
%       correspond to indices \Omega (values are irrevelant).
% Returns
% ----------
% err_fro: double. Frobenius distance ||M-N||_F.
%
% =========================================================================
% Author: Christian Kuemmerle, 2019-2020.
% err_fro = norm_frob_compact(M,sps_plc)^2;

UU=M.U'*N.U;
VV=N.V'*M.V;
m=length(M.res_range);
setSval(sps_plc,M.res_range,m);
M4U=sps_plc'*N.U; % d2 x r
M4V=sps_plc*N.V;  % d1 x r 
setSval(sps_plc,N.res_range,m);
N4U=M.U'*sps_plc; % d2 x r
N4V=sps_plc*M.V;  % d1 x r
M2V=M.Gam2*N.V;
N2V=N.Gam2*M.V;
M3U=M.Gam3'*N.U;
N3U=N.Gam3'*M.U;

lrge_part = norm_frob_compact(M,sps_plc)^2+norm_frob_compact(N,sps_plc)^2 - 2*real(trace(M.Gam1'*UU*N.Gam1*VV));
% sc_prod=trace(M.Gam1'*UU*N.Gam1*VV);
sc_prod=trace((N.Gam1*M2V')*UU)+trace((M.Gam1'*UU)*N2V); %% Potential speed up here? If only the diagonal elements of the outer matrix products are calculated?
sc_prod=sc_prod+trace((M.Gam1'*N3U')*VV)+trace((M3U*N.Gam1)*VV);
sc_prod=sc_prod+trace((N.Gam2*M.Gam2')*UU);
sc_prod=sc_prod+trace(N.Gam2*M4U)+trace(N4U*M.Gam2');
sc_prod=sc_prod+trace(M.Gam3'*N4V)+trace(M4V'*N.Gam3);
sc_prod=sc_prod+trace((M.Gam3'*N.Gam3)*VV);
sc_prod=sc_prod+trace(M.Gam1'*(M.U'*N4V))+trace((M4V'*N.U)*N.Gam1);
sc_prod=sc_prod+trace(M3U*N2V)+trace(M2V'*N3U');
sc_prod=sc_prod+N.res_range*M.res_range';
err_fro = lrge_part - 2*real(sc_prod);
% err_fro = err_fro + norm_frob_compact(N,sps_plc)^2;

err_fro = sqrt(abs(err_fro));

end

