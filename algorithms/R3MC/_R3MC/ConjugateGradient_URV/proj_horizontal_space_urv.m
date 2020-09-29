function[eta] = proj_horizontal_space_urv(model, eta)
    % Project a tangent vector eta onto the horizontal space at model
    
    % Refer "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
    % B. Mishra and R. Sepulchre,
    % Technical report, arXiv:1306.2672, 2013.
    % This implementation is due to
    % Bamdev Mishra <b.mishra@ulg.ac.be>, 2013
    
    
    if ~all(isfield(model,{'RtR','RRt','invRtR','invRRt'}) == 1)
        model.RRt = model.R*model.R';
        model.RtR = model.R'*model.R;
        model.invRRt = eye(size(model.R, 2))/model.RRt;
        model.invRtR = eye(size(model.R, 2))/model.RtR;
    end
    
    
    E = skew((model.U'*eta.U)*model.RRt) + skew(model.R*eta.R');
    F = skew((model.V'*eta.V)*model.RtR) + skew(model.R'*eta.R);
    [Omega1, Omega2] = coupled_lyap(model.R, E, F);% Coupled Lyapunov equations
    
    eta.U = eta.U - (model.U*Omega1);
    eta.R = eta.R - (model.R*Omega2 - Omega1*model.R) ;
    eta.V = eta.V - (model.V*Omega2);
    
    
    %     %% Debug:
    %     norm(skew(model.RRt*(eta.U'*model.U) + eta.R*model.R'), 'fro')
    %     norm(skew(model.RtR*(eta.V'*model.V) - model.R'*eta.R), 'fro')
    
    
    
end

