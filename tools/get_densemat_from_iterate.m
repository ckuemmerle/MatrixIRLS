function X = get_densemat_from_iterate(Xr,sps_plc)
% This functions can be used in two different ways:
% - Given a matrix Xr represented in one of the admissible formats (dense,
% matrix factorization, or compact format of MatrixIRLS), it returns its
% dense representation X (as a (d1 x d2) matrix).
% - Given a cell array of length 'nr_algos' with matrix representations Xr
% in its cells, it calculates a cell array with the respective dense (d1 x
% d2) representation in the cells of a cell array X.
% =========================================================================
% Parameters
% ----------
% Xr:   Option 1: Either of
%       struct. Represents a (d_1 x d_2) matrix X in compact format by 
%               the coefficient matrices
%           Xr.Gam1 = \tilde{Gamma}_1  [(r_c x r_c) matrix]
%           Xr.Gam2 = \tilde{Gamma}_2  [(r_c x d_2) matrix]
%           Xr.Gam3 = \tilde{Gamma}_3  [(d_1 x r_c) matrix]
%           Xr.res_range = r_{k+1}     [(m x 1) vector].
%           Xr.U    = U^{(k)}          [(d_1 x r_c) matrix w/
%                                   orthonormal columns]
%           Xr.V    = V^{(k)}          [(d_2 x r_c) matrix w/
%                                   orthonormal columns]
%       such that  X = U^{(k)}*(\tilde{Gamma}_1*V^{(k)'}+\tilde{Gamma}_2)
%           + \tilde{Gamma}_3*V^{(k)'} + P_{\Omega}^*(r_{k+1}).
%       cell array of size (1 x 2). Represents (d_1 x d_2) matrix X such
%           that X = Xr{1}*Xr{2}'.
%       matrix of size (d_1 x d_2). Dense representation of (d_1 x d_2) 
%           matrix X.
%       Option 2: Cell array of size (1 x nr_algos) with cells as in option
%       1.
%
% sps_plc: (d_1 x d_2) sparse matrix with support of size m. Non-zeros 
%          correspond to indices \Omega (values are irrevelant). This is
%          only relevant if Xr (or any of the cells of Xr) is represented
%          as a struct, see above.
% Returns
% ----------
% X:    (d_1 x d_2) matrix or (1 x nr_algos) cell array with matrices. 
%       Dense representation(s) of the matrix/matrices X.
% =========================================================================
% Author: Christian Kuemmerle, 2020.

if isstruct(Xr) || (iscell(Xr) && size(Xr{1},1) > 1)
    X = get_densemat(Xr,sps_plc);
else
    nr_algos = length(Xr);
    X = cell(1,nr_algos);
    for l=1:nr_algos
        X{l} = get_densemat(Xr{l}{end},sps_plc);
    end
end
end

function Xc = get_densemat(Xrc,sps_plc)
    if isstruct(Xrc)
        Xc = get_densemat_from_compact(Xrc,sps_plc);
    else
        if iscell(Xrc)
            Xc = Xrc{1}*Xrc{2}';
        else
            Xc = Xrc;
        end
    end
end