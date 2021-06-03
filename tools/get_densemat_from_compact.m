function X = get_densemat_from_compact(XX,sps_plc,varargin)
% Given a matrix X represented in the compact format of the Appendix of [1], 
% calculates its dense representation.
% =========================================================================
% Parameters
% ----------
% XX:   struct. Represents a (d_1 x d_2) by the coefficient matrices.
%       XX.Gam1 = \tilde{Gamma}_1  [(r_c x r_c) matrix]
%       XX.Gam2 = \tilde{Gamma}_2  [(r_c x d_2) matrix]
%       XX.Gam3 = \tilde{Gamma}_3  [(d_1 x r_c) matrix]
%       XX.res_range = r_{k+1}     [(m x 1) vector].
%       XX.U    = U^{(k)}          [(d_1 x r_c) matrix w/
%                                   orthonormal columns]
%       XX.V    = V^{(k)}          [(d_2 x r_c) matrix w/
%                                   orthonormal columns],
% sps_plc: (d_1 x d_2) sparse matrix with support of size m. Non-zeros 
%       correspond to indices \Omega (values are irrevelant).
% Returns
% ----------
% X:    (d_1 x d_2) matrix. Dense representation of the matrix such that
%       X = U^{(k)}*(\tilde{Gamma}_1*V^{(k)'}+\tilde{Gamma}_2)
%           + \tilde{Gamma}_3*V^{(k)'} + P_{\Omega}^*(r_{k+1}).
% =========================================================================
% Reference2:
% [1] C. Kuemmerle, C. Mayrink Verdun, "Escaping Saddle Points in 
% Ill-Conditioned Matrix Completion with a Scalable Second Order Method", 
% ICML 2020 Workshop "Beyond First Order Methods in ML Systems".
% [2] C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method for 
% Ill-Conditioned Matrix Completion from Few Samples", ICML 2021.
% =========================================================================
% Author: Christian Kuemmerle, 2019-2020.
if nargin > 3
    ind=varargin{1};
    setSval_ind(sps_plc,ind,XX.res_range,length(XX.res_range));
else
    setSval(sps_plc,XX.res_range,length(XX.res_range));
end
X=(XX.U*XX.Gam1)*XX.V'+XX.U*XX.Gam2+XX.Gam3*XX.V'+sps_plc;

end