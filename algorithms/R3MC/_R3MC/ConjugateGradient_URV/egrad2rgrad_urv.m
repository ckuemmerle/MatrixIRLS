function[rgrad] = egrad2rgrad_urv(model, egrad)
    % egrad is the Euclidean gradient at the point model
    % rgrad is the Riemannian gradient
    
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
    
    
    SSU = model.RRt;
    ASU = 2*symm(SSU*(egrad.R*model.R'));
    SSV = model.RtR;
    ASV = 2*symm(SSV*(egrad.R'*model.R));
    
    [BU, BV] = tangent_space_lyap(model.R, ASU, ASV); % Lyapunov equations on the tangent space
    
    % Riemannian gradient
    rgrad.U = (egrad.U - model.U*BU)*model.invRRt;
    rgrad.V = (egrad.V - model.V*BV)*model.invRtR;
    rgrad.R = egrad.R;
    
    
    % %% Debug: the Riemannian gradient should be on the Horizontal plane
    % norm(skew(model.RRt*(rgrad.U'*model.U) + rgrad.R*model.R'), 'fro')
    % norm(skew(model.RtR*(rgrad.V'*model.V) - model.R'*rgrad.R), 'fro')
    
    
end