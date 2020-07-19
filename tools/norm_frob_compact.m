function val = norm_frob_compact(M,sps_plc,mode)
% Calculates the Frobenius norm ||M||_F of a matrix M given in the the
% compact format
% M_mat =  U^{(k)}*(\tilde{Gamma}_1*V^{(k)'}+\tilde{Gamma}_2)
%          + \tilde{Gamma}_3*V^{(k)'} + P_{\Omega}^*(r_{k+1}).
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
%                                   orthonormal columns
% sps_plc: (d_1 x d_2) sparse matrix with support of size m. Non-zeros 
%       correspond to indices \Omega (values are irrevelant).
% mode: string of characters. If
%        = 'standard': As above.
%        = 'lsqr': (to be updated): M = {M1,M2,M3,M4} such that
%           Mat = [M1         ,  M2*V_perp;
%                  U_perp'*M3 ,  U_perp'*\mathcal{A}'(M4)*V_perp],
%          and U_perp resp. V_perp are orthonormal bases of the
%          spaces that are orthogonal on the columns of U and
%          V, respectively.
% Returns
% ----------
% err_fro: double. Frobenius norm ||M||_F.

%              mode: Indicates if the format of M is such that:
%                    = 'standard': M = {M1,M2,M3,M4} such that
%                                  Mat = U*M1*V' +
%                                  U*M2'+M3*V'+\mathcal{A}'(M4),
%                      and \mathcal{A} is the operator that maps to the
%                      entries represented by sps_plc.
%                    = 'lsqr': M = {M1,M2,M3,M4} such that
%                       Mat = [M1         ,  M2*V_perp;
%                              U_perp'*M3 ,  U_perp'*\mathcal{A}'(M4)*V_perp],
%                      and U_perp resp. V_perp are orthonormal bases of the
%                      spaces that are orthogonal on the columns of U and
%                      V, respectively.
% =========================================================================
% Author: Christian Kuemmerle, 2019-2020.
if nargin == 2
    mode = 'standard';
end
if strcmp(mode,'standard')
    m=length(M.res_range);
    setSval(sps_plc,M.res_range,m);
    M4U=sps_plc'*M.U; % d2 x r
    M4V=sps_plc*M.V;  % d1 x r 
    temp=M.Gam1(:);
    val=sum(temp'*temp);
    val=val+2*real(trace(M.Gam1*(M.V'*M.Gam2')));
    val=val+2*real(trace((M.Gam3'*M.U)*M.Gam1));
    temp=M.Gam2(:);
    val=val+sum(temp'*temp);
    val=val+2*real(trace(M.Gam2*M4U));
    val=val+2*real(trace(M4V'*M.Gam3));
    temp=M.Gam3(:);
    val=val+sum(temp'*temp);
    val=val+2*real(trace(M.Gam1*(M4V'*M.U)));
    val=val+2*real(trace((M.Gam3'*M.U)*(M.Gam2*M.V)));
    temp=M.res_range(:);
    val=val+sum(temp'*temp);
    val=sqrt(val);
elseif strcmp(mode,'lsqr')
    m = length(M.res_range);
    setSval(sps_plc,M.res_range,m);
    UAstM4 = M.U'*sps_plc;
    UAstM4V = UAstM4*M.V;
    AstM4V = sps_plc*M.V; 
    
    temp=M.Gam2(:);
    val2=sum(temp'*temp);
    M2V =M.Gam2*M.V;
    temp=M2V(:);
    val2=val2-sum(temp'*temp);
    
    temp=M.Gam3(:);
    val3=sum(temp'*temp);
    UM3 =M.U'*M.Gam3;
    temp=UM3(:);
    val3=val3-sum(temp'*temp);
    
    temp=M.res_range(:);
    val4=sum(temp'*temp);
    temp=UAstM4(:);
    val4=val4-sum(temp'*temp);
    temp=UAstM4V(:);
    val4=val4+sum(temp'*temp);
    temp=AstM4V(:);
    val4=val4-sum(temp'*temp);
    
    temp=M.Gam1(:);
    val=sum(temp'*temp)+val2+val3+val4;
    val=sqrt(val);
end
end

