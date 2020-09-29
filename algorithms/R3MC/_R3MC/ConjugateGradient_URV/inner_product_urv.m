function[ip] = inner_product_urv(model, eta, xi)
    % Inner product between tangent vectors eta and xi at the point model
    
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
    
    ip = trace(model.RRt*(eta.U'*xi.U)) + trace(model.RtR*(eta.V'*xi.V)) + trace(eta.R'*xi.R);
end
