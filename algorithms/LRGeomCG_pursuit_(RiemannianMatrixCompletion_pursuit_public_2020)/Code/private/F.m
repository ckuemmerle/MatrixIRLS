function [f, df] = F(prob,x, dx)
%F Cost function 
%
% df is the gradient of f at x, since x 

if nargout==2
    if nargin<3
        error('need dx')
    end
else
    dx = 0;
end


e = x.on_omega - prob.data;
f = 0.5*(e'*e);

if prob.mu>0
  f = f + prob.mu*0.5*norm(1./x.sigma)^2;
end




if nargout==2
    updateSval(prob.temp_omega, e, length(e));
    prob.temp_omega = prob.temp_omega*1;
    T_v = prob.temp_omega*x.V;
    T_u = prob.temp_omega*x.U;
    df = trace( dx.M'*(x.U'*T_v) + dx.Up'*T_v + T_u'*dx.Vp);     
end