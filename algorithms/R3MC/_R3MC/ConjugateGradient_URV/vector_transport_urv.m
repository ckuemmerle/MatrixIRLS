function[eta_new] = vector_transport_urv(model, eta, model_new)
    % Transport vector eta at model to model_new
    
    % Refer "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
    % B. Mishra and R. Sepulchre,
    % Technical report, arXiv:1306.2672, 2013.
    % This implementation is due to
    % Bamdev Mishra <b.mishra@ulg.ac.be>, 2013
    
    
    if ~all(isfield(model_new,{'RtR','RRt','invRtR','invRRt'}) == 1)
        model_new.RRt = model_new.R*model_new.R';
        model_new.RtR = model_new.R'*model_new.R;
        model_new.invRRt = eye(size(model_new.R, 2))/model_new.RRt;
        model_new.invRtR = eye(size(model_new.R, 2))/model_new.RtR;
    end
    
    
    
    ip_old = inner_product_urv(model, eta, eta);
    
    
    % Project eta on the tangent sapce of model_new
    eta = proj_tangent_space_urv(model_new, eta);
    % %% Debug:
    %
    % norm(symm(model_new.U'*eta.U), 'fro')
    % norm(symm(model_new.V'*eta.V), 'fro')
    
    
    
    
    % Project eta on to the horizontal space
    eta_new = proj_horizontal_space_urv(model_new, eta);
    % %% Debug:
    %
    % norm(skew(model_new.RRt*(eta_new.U'*model_new.U) + eta_new.R*model_new.R'), 'fro')
    % norm(skew(model_new.RtR*(eta_new.V'*model_new.V) - model_new.R'*eta_new.R), 'fro')
    
    
    ip_new = inner_product_urv(model_new, eta_new, eta_new);
    
    scaling = (ip_old / ip_new);
    
    eta_new.U = scaling * eta_new.U;
    eta_new.R = scaling * eta_new.R;
    eta_new.V = scaling * eta_new.V;
    
    
end
