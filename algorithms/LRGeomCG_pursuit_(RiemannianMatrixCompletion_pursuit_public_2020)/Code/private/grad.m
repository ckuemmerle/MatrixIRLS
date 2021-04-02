function g = grad(prob,x)
%GRAD   Computes the gradient on the manifold
%
%  Computes the gradient at a point x of the cost function on 
%  the manifold.



d = x.on_omega - prob.data;
updateSval(prob.temp_omega, d, length(d));
prob.temp_omega = prob.temp_omega*1;

T = prob.temp_omega*x.V;
g.M = x.U'*T; 
g.Up = T - x.U*(x.U'*T);
g.Vp = prob.temp_omega'*x.U; 
g.Vp = g.Vp - x.V*(x.V'*g.Vp);

    
if prob.mu>0
  g.M = g.M - prob.mu*diag(x.sigma.^(-3));
end

