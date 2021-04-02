function diff_norm = compute_iterate_diff_norm(Xc,Xc_,sps_plc,mode)
%compute_iterate_diff_norm This function computes the Frobenius norm of the
%difference of matrices represented in the "low-rank + sparse" format of
%the paper
% [1] C. Kuemmerle, C. Mayrink Verdun, "A Scalable Second Order Method 
% for Ill-Conditioned Matrix Completion from Few Samples".

if strcmp(mode,'exact') % this is exact (does not suffer from numerical errors 
    % very much), but too slow for large problems as full (d1 x d2)
    % matrices are computed
    Xcfull = get_densemat_from_compact(Xc,sps_plc);
    Xc_full = get_densemat_from_compact(Xc_,sps_plc);
    diff_norm = norm(Xcfull-Xc_full,'fro');
elseif strcmp(mode,'efficient') % this is computationally efficient as it 
    % does not compute full matrices, but as of Apr 2021, its
    % implementation suffers from issues in the numerical precision.
    diff_norm = get_frob_diff_compact(Xc,Xc_,sps_plc);
elseif strcmp(mode,'proxy') % this computes a proxy that needs only O(m)
    % operations and is of the same magnitude as the exact value.
    % Sufficient for most purposes (such as stopping conditions).
    diff_norm = 5*norm(Xc.res_range-Xc_.res_range);
else
    error('Error in the computation of the norm of differences of iterates.')
end
end

