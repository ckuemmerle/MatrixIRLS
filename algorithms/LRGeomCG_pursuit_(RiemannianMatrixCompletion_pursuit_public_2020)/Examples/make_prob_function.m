function prob = make_prob_function(fun_f, bnds, n1,n2, r,OS,testing,opt_S_fast)
% Make a matrix completion problem X= L*S*R' with random L and R and
% specified S.
% 
% Arguments: 
%   [n1,n2] = size of X
%   r = rank of X        
%   s = diag of S
%   OS = oversampling ratio in Omega wrt dof in X
%   testing = true if testing set (same size as Omega)

if nargin==7
    opt_S_fast = false;
end

prob.n1 = n1;
prob.n2 = n2;
prob.r = r;
if prob.r >= min(prob.n1, prob.n2)
    error('The rank of the matrices in the reconstruction manifold should be stricly smaller than number of rows and columns.')
end

x = linspace(bnds(1), bnds(2), n1);
y = linspace(bnds(3), bnds(4), n2);

dof = r*(n1+n2-r);

% make random sampling, problem and initial guess
samples = floor(OS * dof);
Omega = make_rand_Omega(n1,n2,samples); % training set

prob.Omega = Omega;
[prob.Omega_i, prob.Omega_j] = ind2sub([prob.n1,prob.n2], Omega);
prob.m = length(Omega);

prob.data = fun_f(x(prob.Omega_i), y(prob.Omega_j))';
prob.norm_M_Omega = norm(prob.data);

% [X,Y] = meshgrid(x,y);
% Z = fun_f(X,Y);
% d2 = partXY(Z', eye(size(Z)),prob.Omega_i, prob.Omega_j,prob.m)'; 
% norm(d2 - prob.data)
% pause


prob.temp_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data*1,prob.n1,prob.n2,prob.m);
  
% regularization parameter
%prob.mu = 1e-14;
prob.mu = 0;


if testing
    Gamma = make_rand_Omega(n1,n2,samples); % testing set

    prob.has_testing = true;
    prob.Gamma = Gamma;
    [prob.Gamma_i, prob.Gamma_j] = ind2sub([prob.n1,prob.n2], Gamma);
    prob.test_data = fun_f(x(prob.Gamma_i), y(prob.Gamma_j))';    
    %prob.test_data = partXY(Z', eye(size(Z)),prob.Gamma_i, prob.Gamma_j,prob.m)'; 
else
    prob.has_testing = false;
end

if opt_S_fast
    inds =  rand_indices(length(Omega), 2*(r)^2 + 5);
    
    Omega_S = Omega( inds );
    %Omega = make_rand_Omega(n1,n2,3*(r+5)^2); % a little bigger than dof in S

    prob.opt_S_fast = true;
    prob.opt_S_fast.Omega = Omega_S;
    [prob.opt_S_fast.Omega_i, prob.opt_S_fast.Omega_j] = ind2sub([prob.n1,prob.n2], Omega_S);
    prob.opt_S_fast.m = length(Omega_S);
    prob.opt_S_fast.data = fun_f(x(prob.opt_S_fast.Omega_i), y(prob.opt_S_fast.Omega_j))';
    
    prob.opt_S_fast.temp_omega = sparse(prob.opt_S_fast.Omega_i, prob.opt_S_fast.Omega_j, ...
            prob.opt_S_fast.data*1,prob.n1,prob.n2,prob.opt_S_fast.m);    
else
    prob.opt_S_fast = false;
end
end

function omega = rand_indices(n,m)
% random uniformly distributed on [0,n]
omega = ceil(rand(m, 1) * n);
omega = unique(omega);
while length(omega) < m    
    omega = [omega; ceil(rand(m-length(omega), 1)*n);];
    omega = unique(omega);
end
omega = omega(1:m);
omega = sort(omega);

end
    
