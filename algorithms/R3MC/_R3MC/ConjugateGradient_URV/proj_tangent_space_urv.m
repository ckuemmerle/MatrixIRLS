function eta = proj_tangent_space_urv(model, eta)
    % Projecting a vector eta onto the tangent space of URV
    %
    
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
    
    
    SSU = model.RRt; %model.RRt;
    ASU = 2*symm(SSU*(model.U'*eta.U)*SSU);
    SSV = model.RtR; %model.RtR;
    ASV = 2*symm(SSV*(model.V'*eta.V)*SSV);
    
    [BU, BV] = tangent_space_lyap(model.R, ASU, ASV); % Lyapunov equations on the tangent space
    
    % Projection
    eta.U = eta.U - model.U*(BU*model.invRRt);
    eta.V = eta.V - model.V*(BV*model.invRtR);
    
    
    % %% Debug:
    %
    % norm(symm(model.U'*eta.U), 'fro')
    % norm(symm(model.V'*eta.V), 'fro')
    
    
end