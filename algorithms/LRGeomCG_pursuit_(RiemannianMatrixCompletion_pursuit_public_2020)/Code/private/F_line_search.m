function [f, xt, store, grad_t, df] = F_line_search(prob,x,dir, t, store)
%Evaluate the objective function along a retracted ray
%   f = f_x( t ) := F( R_x( t*dir ) )
% where dir is tangent vector at x.
%
% Returns:
%     f = f_x( t )
%     xt = R_x( t*dir )
%     grad_t is the Riemannian gradient of F at F( R_x ( t*dir ) )
%     df is gradient of f_x at t. 
%     store: computationa results specific to prob, x, dir. Can be reused
%            for different t
%



k = size(x.V,2);
Zero = zeros(k);

if nargin==4 || isempty(store)
    % qr of Vp and Up only
    [store.U_t,store.Ru] = qr(dir.Up,0);
    [store.V_t,store.Rv] = qr(dir.Vp,0);
end

% new point xt
Ru_M_Rv = [x.S+t*dir.M t*store.Rv'; t*store.Ru Zero];

    [UU,SS,VV] = svd(Ru_M_Rv);


xt.U = [x.U store.U_t]*UU(:,1:k);
xt.V = [x.V store.V_t]*VV(:,1:k);
xt.S = SS(1:k,1:k);

xt = prepx(prob, xt);

% Function value at xt
e = xt.on_omega - prob.data;
f = 0.5*(e'*e);

if prob.mu>0
  f = f + prob.mu*0.5*norm(1./diag(xt.S))^2;
end

if nargout>3
% Riemannian gradient at xt
updateSval(prob.temp_omega, e, length(e));
prob.temp_omega = prob.temp_omega*1;

T_v = prob.temp_omega*xt.V;
T_u = prob.temp_omega'*xt.U;

grad_t.M = xt.U'*T_v; 
grad_t.Up = T_v - xt.U*grad_t.M;
grad_t.Vp = T_u - xt.V*grad_t.M';

if prob.mu>0
  grad_t.M = grad_t.M - prob.mu*diag(diag(xt.S).^(-3));
end
end


if nargout>4
% derivative of f_x
SS2 = SS(k+1:2*k,k+1:2*k);
SS1 = SS(1:k,1:k);

D = UU'*[dir.M store.Rv'; store.Ru Zero]*VV;
D11 = D(1:k,1:k);
D12 = D(1:k,k+1:2*k);
D21 = D(k+1:2*k,1:k); 

ss1 = diag(SS1).^2;
ss2 = diag(SS2).^2;
H12 = (D12*SS2+SS1*D21')./(ones(k,1)*ss2'-ss1*ones(1,k));
K12 = SS1\(H12*SS2-D12);

H11 = (D11-D11')./(ones(k,1)*ss1'+ss1*ones(1,k));
K11 = Zero; 

dS1 = D11 + SS1*K11 - H11*SS1;

dU1 = UU*[H11; -H12'];
dV1 = VV*[K11; -K12'];

dx_t.M = dS1;
dx_t.Up = [x.U store.U_t]*dU1*SS(1:k,1:k);
dx_t.Vp = [x.V store.V_t]*dV1*SS(1:k,1:k)';


df = trace( dx_t.M'*(xt.U'*T_v) + dx_t.Up'*T_v + dx_t.Vp'*T_u);     
end