function[BU, BV] = tangent_space_lyap(R, E, F)
    % We intent to solve     RR^T  BU + BU RR^T  = E
    %                        R^T R BV + BV R^T R = F
    %
    % This can be solved using two calls to the Matlab lyap.
    % However, we can still have a more efficient implementations as shown
    % below.
    
    % Refer "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
    % B. Mishra and R. Sepulchre,
    % Technical report, arXiv:1306.2672, 2013.
    % This implementation is due to
    % Bamdev Mishra <b.mishra@ulg.ac.be>, 2013
    
    
    [U, Sigma, V] = svd(R);
    E_mod = U'*E*U;
    F_mod = V'*F*V;
    b1 = E_mod(:);
    b2 = F_mod(:);
    
    
    r = size(Sigma, 1);
    sig = diag(Sigma); % all the singular values in a vector
    sig1 = sig*ones(1, r); % columns repeat
    sig1t = sig1'; % rows repeat
    s1 = sig1(:);
    s2 = sig1t(:);
    
    % The block elements
    a = s1.^2 + s2.^2; % a column vector
    
    % solve the linear system of equations
    cu = b1./a; %a.\b1;
    cv = b2./a; %a.\b2;
    
    % devectorize
    CU = reshape(cu, r, r);
    CV = reshape(cv, r, r);
    
    % Do the similarity transforms
    BU = U*CU*U';
    BV = V*CV*V';
    
    % %% debug
    %
    % norm(R*R'*BU + BU*R*R' - E, 'fro');
    % norm((Sigma.^2)*CU + CU*(Sigma.^2) - E_mod, 'fro');
    % norm(a.*cu - b1, 'fro');
    %
    % norm(R'*R*BV + BV*R'*R - F, 'fro');
    %
    % BU1 = lyap(R*R', - E);
    % norm(R*R'*BU1 + BU1*R*R' - E, 'fro');
    %
    % BV1 = lyap(R'*R, - F);
    % norm(R'*R*BV1 + BV1*R'*R - F, 'fro');
    %
    % % as accurate as the lyap
    % norm(BU - BU1, 'fro')
    % norm(BV - BV1, 'fro')
    
end