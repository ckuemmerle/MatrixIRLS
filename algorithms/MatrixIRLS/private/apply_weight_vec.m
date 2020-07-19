function gam_out = apply_weight_vec(gam,d1,d2,weight_vec,increase_antisymmetricweights,varargin)
%apply_weight_vec This function just applies an elementwise multiplication
% of gam with the elements of weight_vec in the case of most 'type_mean'.
% However, for cases where weights are chosen according to a symmetric-skew-
% symmetric splitting as type_mean = 'Hessianupper_2', the multiplication
% is not diagaonal, but takes this splitting into account.
if increase_antisymmetricweights
    r=round((-(d1+d2)+sqrt((d1+d2)^2+4*length(gam)))./2);
    if not(length(gam) == r*(d1+d2+r))
        error('Error in the dimensionality.')
    end
    
    gam_Tkanti = get_antisymmTkbasis_from_Tk(gam,d1,d2,r);
    gam_Tkanti = weight_vec.*gam_Tkanti;
    gam_weighted = get_Tkbasis_from_antisymmTk(gam_Tkanti,d1,d2,r);
    gam_out=gam_weighted;
else
        gam_out = weight_vec.*gam;
end
end

