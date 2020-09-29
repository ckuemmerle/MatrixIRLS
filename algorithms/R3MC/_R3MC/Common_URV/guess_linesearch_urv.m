function[xmin] = guess_linesearch_urv(model, errors, searchDir, data_ls)
    % Full linearized stepsize search in the direction searchDir at the point model
    
    % Refer "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
    % B. Mishra and R. Sepulchre,
    % Technical report, arXiv:1306.2672, 2013.
    % This implementation is due to
    % Bamdev Mishra <b.mishra@ulg.ac.be>, 2013
    
    
    n = data_ls.nentries;
    searchDir.U = n*searchDir.U;
    searchDir.V = n*searchDir.V;
    searchDir.R = n*searchDir.R;
    
    term1 = (searchDir.U*model.R + model.U*searchDir.R);
    term2 = (searchDir.U*searchDir.R);
    
    s_star = partXY([model.U*model.R  term1]', [searchDir.V  model.V]', data_ls.rows, data_ls.cols, n)';
    s_star_1 = partXY([term1 term2]', [searchDir.V  model.V ]', data_ls.rows, data_ls.cols, n)';
    s_star_urv = partXY(term2', searchDir.V', data_ls.rows, data_ls.cols, n)';
    
    
    
    % A degree 6 polynomial
    a_6 = (s_star_urv'*s_star_urv); % sum(s_star_urv.^2);
    a_5 = 2*(s_star_urv'*s_star_1);
    a_4 = 2*(s_star'*s_star_urv) + (s_star_1'*s_star_1); % 2*sum(s_star .* s_star_urv) + sum(s_star_1.^2);
    a_3 = 2*(s_star'*s_star_1) + 2*(errors'*s_star_urv); % 2*sum(s_star .* s_star_1) + 2*sum(errors .* s_star_urv);
    a_2 = 2*(s_star_1'*errors) + (s_star'*s_star);% 2*sum(s_star_1 .* errors) + sum(s_star .^2);
    a_1 = 2*(s_star'*errors); % 2*sum(s_star .* errors);
    a_0 = (errors'*errors);% sum(errors.^2);
    % old_cost = a_0 / n; % old cost
    
    coefficients_poly = [a_6 a_5 a_4 a_3 a_2 a_1 a_0];
    
    % Solving the derivative of the degree 6 polynomial to find the global
    % minimum
    coefficients = [6*a_6 , 5*a_5, 4*a_4, 3*a_3, 2*a_2, a_1];
    xs = roots(coefficients);
    
    real_roots = xs(imag(xs)==0);
    real_roots = 0.5*(real_roots + abs(real_roots)); % Positive roots
    
    fun_val = polyval(coefficients_poly, real_roots);
    [~, idx] = min(fun_val);
    
    xmin = real_roots(idx); % Global minimum
    xmin = n*xmin; % Rescaling
end






