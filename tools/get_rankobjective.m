function objective_val = get_rankobjective(sing,eps,p,mode)
%get_rankobjective This function evaluates the logdet or Schatten-p
%objective associated to a given vector of (partial) singular values and
%given a value of the epsilon smoothing eps.
%  Input :               
%                        sing = vector of singular values with values
%                               larger than the smoothing parameter epsilon
%                         eps = value of smoothing parameter epsilon
%                           d = minimum of matrix dimensions d1 and d2
%                           p = quasi-norm parameter: Schatten-p objective
%                               if p > 0, logdet objective for p = 0.
%                        mode = 'objective_thesis': Corresponds to using objectives
%                           of the type:
%                           log(max(sigma,eps))+0.5.*(min(sigma,eps)^2/eps^2-1).
%                           (as used in the Ph.D. thesis [Kuemmerle
%                           2019])
%                             = 'pluseps': Corresponds to using objectives
%                           of the type:
%                           log(max(sigma,eps)+eps)+0.25.*(min(sigma,eps)^2/eps^2-1)
%                             = 'pluseps_squared':Corresponds to using objectives
%                           of the type:
%                           log(max(sigma,eps)^2+eps^2)+0.5.*(min(sigma,eps)^2/eps^2-1)
r_c=length(sing);
switch mode
    case 'objective_thesis'
        if p == 0
            objective_val = sum(log(max(sing,eps))+0.5.*(min(sing,eps).^2./eps.^2-1));
        else
            objective_val = sum(max(sing,eps).^p+(p/2).*(min(sing,eps).^2./eps.^(2-p)-1));
        end
    case 'pluseps'
        if p == 0
            objective_val = sum(log(max(sing,eps)+eps)+0.25.*(min(sing,eps).^2./eps.^2-1));
        else
            objective_val = sum((max(sing,eps)+eps).^p+(p*(2*eps)^p).*(min(sing,eps).^2./(4*eps.^2)-0.25));
        end
    case 'pluseps_squared_max'
        if p == 0
            objective_val = sum(log(max(sing,eps).^2+eps^2)+0.5.*(min(sing,eps).^2./eps.^2-1));
        else
            objective_val = sum(2.*(max(sing,eps).^2+eps^2).^(p/2)+(p*(2*eps^2)^(p/2)).*(min(sing,eps).^2./(2*eps.^2)-0.5));
        end
    case 'pluseps_squared'
        if p == 0
            objective_val = sum(log(sing.^2+eps^2));
        else
            objective_val = sum(2.*(sing.^2+eps^2).^(p/2));
        end
end

end

