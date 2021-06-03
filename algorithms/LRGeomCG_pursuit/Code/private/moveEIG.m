function [z, dz] = moveEIG(prob, x,h,t)
%MOVEEIG  Retract a point x+t*h (minimal distance)
%
%  z = moveEIG(sys, x,h,t) retracts x+t*h to the manifold where
%   sys is a system,
%     x is a point of sys,
%     h is a tangent vector of TxM,
%     t is a distance (scalar).
%  Retraction is the shortest distance to the manifold.
%
%  Implemented as a partitioned eigenvalue decomposition with 
%  complexity O(k^3 n).


k = size(x.V,2);
Zero = zeros(k);

% qr of Vp and Up only
[U_t,Ru] = qr(h.Up,0);
[V_t,Rv] = qr(h.Vp,0);


Ru_M_Rv = [diag(x.sigma)+t*h.M t*Rv'; t*Ru Zero];
[UU,SS,VV] = svd(Ru_M_Rv);

z.U = [x.U U_t]*UU(:,1:k);
z.V = [x.V V_t]*VV(:,1:k);
z.sigma = diag(SS(1:k,1:k));

z = prepx(prob, z);

if nargout==2
   SS2 = SS(k+1:2*k,k+1:2*k);
   SS1 = SS(1:k,1:k);
    
   D = UU'*[h.M Rv'; Ru Zero]*VV;
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
   
   dz.M = dS1;
   dz.Up = [x.U U_t]*dU1*SS(1:k,1:k);
   dz.Vp = [x.V V_t]*dV1*SS(1:k,1:k)';
end

end
