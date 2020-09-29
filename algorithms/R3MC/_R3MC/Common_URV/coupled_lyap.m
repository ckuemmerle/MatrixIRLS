function[Omega1, Omega2] = coupled_lyap(R, E, F)
    % We intent to solve the coupled system of Lyapunov equations
    %
    % RR^T Omega1 + Omega1 RR^T  - R Omega2 R^T = E
    % R^T R Omega2 + Omega1 R^T R  - R^T Omega2 R = F
    %
    % Below is an efficient implementation
    
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
    sig = diag(Sigma); % All the singular values in a vector
    sig1 = sig*ones(1, r); % Columns repeat
    sig1t = sig1'; % Rows repeat
    s1 = sig1(:);
    s2 = sig1t(:);
    
    % The block elements
    a =  s1.^2 + s2.^2; % a column vector
    c = s1.*s2;
    
    
    % Solve directly using the formula
    % A = diag(a);
    % C = diag(c);
    % Ay1 - Cy1 = b1
    % Ay2 - Cy2 = b2;
    % y1_sol = (A*(C\A) - C) \ (b2 + A*(C\b1));
    % y2_sol = A\(b2 + C*y1_sol);
    
    y1_sol = (b2 + (a.*b1)./c) ./ ((a.^2)./c - c);
    y2_sol = (b2 + c.*y1_sol)./a;
    
    
    % devectorize
    Omega1 = reshape(y1_sol, r, r);
    Omega2 = reshape(y2_sol, r, r);
    
    
    % Do the similarity transforms
    Omega1 = U*Omega1*U';
    Omega2 = V*Omega2*V';
    
    
    % %% debug whether we have the right solution
    % norm(R*R'*Omega1 + Omega1*R*R'  - R*Omega2*R' - E, 'fro')
    % norm(R'*R*Omega2 + Omega2*R'*R  - R'*Omega1*R - F, 'fro')
    
    
    
end