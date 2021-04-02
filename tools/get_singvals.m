function sings = get_singvals(X,sps_plc,r)
% This function calculates r singular values of matrix X.
% =========================================================================
% Parameters
% ----------
% Returns
% ----------
% =========================================================================
% Author: Christian Kuemmerle, 2020.
    if isstruct(X)
        d = min(size(X.U,1),size(X.V,1));
        X_handle = get_handle_from_compact(X,sps_plc);
        [~, sing_c, ~]=bksvd_mod(X_handle,...
                min(d,r), 20);
        sings=diag(sing_c);
    else
        if iscell(X)
            d1 = size(X{1},1);
            d2 = size(X{2},1);
            d = min(d1,d2);
            Xc_hdl = @(x) X{1}*(X{2}'*x);
            Xct_hdl = @(y) X{2}*(X{1}'*y);
            X_handle = {Xc_hdl,Xct_hdl,d1,d2};
            [~, sing_c, ~]=bksvd_mod(X_handle,...
                min(d,r), 20);
            sings=diag(sing_c);
        else
            sings=svds(X,r);
        end
    end


end

function  X_c_handle = get_handle_from_compact(X,sps_plc)
d1 = size(X.U,1);
d2 = size(X.V,1);
setSval(sps_plc,X.res_range,length(X.res_range));
Xc_hdl  = @(x) (X.U*(X.Gam1*(X.V'*x) +(X.Gam2*x))...
    + X.Gam3*(X.V'*x) +sps_plc*x);
Xct_hdl = @(y) (X.V*(X.Gam1'*(X.U'*y)+X.Gam3'*y)...
    + X.Gam2'*(X.U'*y) +sps_plc'*y);
X_c_handle={Xc_hdl,Xct_hdl,d1,d2};
end