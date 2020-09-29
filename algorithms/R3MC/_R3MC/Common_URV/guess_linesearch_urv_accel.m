function[xmin] = guess_linesearch_urv_accel(model, errors, searchDir, data_ls)
    % Find the approximated linearized stepsize in the direction searchDir at the point
    % model.
    
    % Refer "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
    % B. Mishra and R. Sepulchre,
    % Technical report, arXiv:1306.2672, 2013.
    % This implementation is due to
    % Bamdev Mishra <b.mishra@ulg.ac.be>, 2013
    
    n = data_ls.nentries;
    
    % Accelerating by reformulating as one factor factorization which reduces to
    % degree 2 from degree 6
    
    term1 = (searchDir.U*model.R + model.U*searchDir.R);
    s_star = partXY([model.U*model.R  term1]', [searchDir.V  model.V]', data_ls.rows, data_ls.cols, n)';
    
    term_num = (errors'*s_star);
    term_den = (s_star'*s_star);
    
    xmin = - term_num/term_den;
    xmin = max(0, xmin);
    
end

