function x = prepx(prob, x)
%PREPX  Prepare a point
%
%  Prepare a point x on the manifold to include extra information by
%  computing
%   

if isempty(x.S)
    x.on_omega = zeros(size(prob.Omega));
else
    x.on_omega = partXY(x.S'*x.U',x.V', prob.Omega_i, prob.Omega_j, prob.m)'; 
end


