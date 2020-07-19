function [Sc,Sc_eps] = set_weight_infovec(sing,eps,p,objective,d,steep_weights)
%set_weight_infovec This function calculates the values from which the weight matrices
% are defined for IRLS algorithms, given a choice of an objective function 
% corresponding to variable 'objective'.
if nargin == 4
    d = [];
    steep_weights = 0;              
end
if steep_weights
    eps_steep=eps./d^p;
else
    eps_steep=eps;
end

switch objective
    case 'objective_thesis'
        Sc = max(sing,eps).^(2-p);
        Sc_eps    = eps_steep^(2-p); 
    case 'pluseps'
        Sc = (max(sing,eps)+eps).^(1-p).*max(sing,eps);
        Sc_eps    = 2^(1-p)*eps_steep^(2-p);
    case 'pluseps_squared_max'
        Sc = ((max(sing,eps).^2+eps^2).^(1-p/2))./2;
        Sc_eps    = 2^(-p/2)*eps_steep^(2-p);    
    case 'pluseps_squared'
        Sc = (sing.^2+eps^2).^(1-p/2)./2;
        Sc_eps    = eps_steep.^(2-p)/4;    
    otherwise
        error('Something is wrong with the choice of the objective.')
end
        
end

