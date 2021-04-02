function [Xt,Y_transp] = transpVect_differentiated(prob, x,h,d, type)
%function [Xt, Y_transp] = transpVect_differentiated(x,z, y)

%% z is the foot of ht
if nargin == 5 && type==1
    error('not implemented')    
end

z = moveEIG(prob,x,d,1);


% vector y belonging to x, transported to R_x(z)

ip_old = ip(x,h,h);

[n,r] = size(d.Up);

[Up,R] = qr([d.Up h.Up],0); Rz = R(:,1:r); Ry = R(:,r+1:2*r);
[Vp,S] = qr([d.Vp h.Vp],0); Sz = S(:,1:r); Sy = S(:,r+1:2*r);

% X = x.U*x.S*x.V';
% Z = x.U*d.M*x.V' + d.Up*x.V' + x.U*d.Vp';
% Y = x.U*h.M*x.V' + h.Up*x.V' + x.U*h.Vp';

% t = 2;
% T = X+Z+t*Y;
% TT = [x.U Up] * [x.S+d.M+t*h.M Sz'+t*Sy'; Rz+t*Ry zeros(2*r)] * [x.V Vp]';
% norm(T-TT)


[n,k] = size(d.Up);
Zero = zeros(2*k);

[UU,SS,VV] = svd([diag(x.sigma)+d.M Sz'; Rz Zero]);
%Xt = [X.U Qu]*UU(:,1:k)*SS(1:k,1:k)*VV(:,1:k)'*[X.V Qv]';
Xt.U = [x.U Up]*UU(:,1:k);
Xt.sigma = diag(SS(1:k,1:k));
Xt.V = [x.V Vp]*VV(:,1:k);



    
    SS2 = SS(k+1:3*k,k+1:3*k);
    SS1 = SS(1:k,1:k);
    
    D = UU'*[h.M Sy'; Ry Zero]*VV;    
    D11 = D(1:k,1:k);
    D12 = D(1:k,k+1:3*k);
    D21 = D(k+1:3*k,1:k);
    D22 = D(k+1:3*k,k+1:3*k);    
    % long version via Sylvester equation, may be more stable?
     Z = lyap(-[zeros(k) SS1; SS1 zeros(k)], SS2, -[D12; D21']);          
     H12 = Z(1:k,1:2*k);
     K12 = Z(k+1:2*k,1:2*k);
%     %H12 = lyap(-SS1*SS1, SS2*SS2, -D12*SS2-SS1*D21'); % is actually a Cauchy matrix    
%     ss1 = diag(SS1).^2;
%     ss2 = diag(SS2).^2;
%     H12 = (D12*SS2+SS1*D21')./(ones(k,1)*ss2'-ss1*ones(1,k));
%    K12 = SS1\(H12*SS2-D12);
    
    H11 = lyap(SS1, -(D11-D11')); % is actually a Cauchy matrix    
%    H11 = (D11-D11')./(ones(k,1)*ss1'+ss1*ones(1,k));
    %H22 = lyap(SS2, -(D22-D22'));
    K11 = zeros(k); 
    %K22 = Zero;
    
    %dS = D + blkdiag(SS1,SS2)*[K11 K12; -K12' K22] - [H11 H12; -H12' H22]*blkdiag(SS1,SS2);
    %dS1 = dS(1:k,1:k); 
    %dS2 = dS(k+1:2*k,k+1:2*k);
    dS1 = D11 + SS1*K11 - H11*SS1;
    
    dU1 = UU*[H11; -H12'];
    dV1 = VV*[K11; -K12'];
    
    %dA = dU1*SS(1:k,1:k)*VV(:,1:k)' + ...
    %    UU(:,1:k)*dS1*VV(:,1:k)' + ...
    %    UU(:,1:k)*SS(1:k,1:k)*dV1';
    %
    %dXt = [X.U Qu]*dA*[X.V Qv]';
    
    Y_transp.M = dS1;
    Y_transp.Up = [x.U Up]*dU1*SS(1:k,1:k);
    Y_transp.Vp = [x.V Vp]*dV1*SS(1:k,1:k)';
    
ip_new = ip(z,Y_transp,Y_transp);
Y_transp = scaleTxM(Y_transp, ip_old/ip_new);

