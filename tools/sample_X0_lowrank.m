function [U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complex,varargin)
% Constructs a random rank-r matrix X0 of dimension (d1 x d2). 
% The randomness comes from Gaussian matrices, parametrized such that
% X0 = U0 * V0' where U0 and V0 have r columns.
% =========================================================================
% Parameters
% ----------
% d1:       number of rows of low-rank matrix X0.
% d2:       number of columns of low-rank matrix X0.
% r:        rank of X0.
% modeX0:   char string or float/int. Indicates the random model
%           for X0.
%           = 1: random model as in the paper [1] (Section 5, X0 = U*S*V' 
%           where U,V are rank-r with i.i.d. standard Gaussian entries
%           and S is (r x r), diagonal and has i.i.d. standard Gaussian 
%           entries.
%           = 2: random model such that X0 = U*V' where U,V are
%           rank-r with i.i.d. standard Gaussian entries.
%           = 'condition_control_linear': X0 = U*S*V' with random U, V with
%           orthonormal columns, and S such that S(1) = cond_nr, S(r) = 1,
%           with _linear_ interpolation in between.
%           = 'condition_control_1/x2': X0 = U*S*V' with random U, V with
%           orthonormal columns, and S such that S(1) = cond_nr, S(r) = 1,
%           with interpolation of the order (1/x^2) in between.
%           = 'condition_control_log': X0 = U*S*V' with random U, V with
%           orthonormal columns, and S such that S(1) = cond_nr, S(r) = 1,
%           with _exponential_ interpolation in between.
% complex:  boolean. If =1, complex matrix X0. If =0, real matrix X0.
% Optional argument: int. Indicates the condition number 'cond_nr', for the
%           random models that involve a control of the condition number.
% Returns
% ----------
% U0:       First factor matrix of X0.
% V0:       Second factor matrix of X0.
% =========================================================================
% References:
% [1] C. Kuemmerle, J. Sigl, "Harmonic Mean Iteratively Reweighted Least 
% Squares for Low-Rank Matrix Recovery", Journal of Machine Learning 
% Research, 19(47):1?49, 2018.
% [2] C. Kuemmerle, C. M. Verdun, "Escaping Saddle Points in 
% Ill-Conditioned Matrix Completion with a Scalable Second Order Method", 
% ICML 2020 Workshop "Beyond First Order Methods in ML Systems".
% =========================================================================
if isfloat(modeX0)
    if (modeX0 == 1 || modeX0 == 2)
        if complex == 1
            U = (1/sqrt(2)).*(randn(d1,r)+1i.*randn(d1,r));
            V = (1/sqrt(2)).*(randn(d2,r)+1i.*randn(d2,r));
            if modeX0 == 1
                S = (1/sqrt(2)).*(randn(1,r)+1i.*randn(1,r));
            elseif modeX0 == 2
                S = ones(1,r);
            end
        else
            U = randn(d1,r);
            V = randn(d2,r);
            if modeX0 == 1
                S = randn(1,r);
            elseif modeX0 == 2
                S = ones(1,r);
            end
        end
    else
        error('This value of "modeX0" is not supported, please choose "1" or "2" or a character string.')
    end
    U0 = U.*S;
elseif strcmp(modeX0,'condition_control_1/x2') || strcmp(modeX0,'condition_control_linear') ...
        || strcmp(modeX0,'condition_control_log') || strcmp(modeX0,'condition_control_log_plateau') ...
        || strcmp(modeX0,'condition_control_log_plateau_end') || strcmp(modeX0,'condition_control_log_plateau_NEW')
    if complex == 1
        error('Complex not yet implemented for matrices with controlled condition number.')
    end
    cond_nr = varargin{1};
%     [U,~,V]=svds(rand(d1,d2),r);
    U = randn(d1,r);
    [U,~,~]=svd(U,'econ');
    V = randn(d2,r);
    [V,~,~]=svd(V,'econ');
%     U = U(:,1:r);
%     V = V(:,1:r);
    if strcmp(modeX0,'condition_control_1/x2')
        fct = @(l) (cond_nr - (1-cond_nr/r^2)/(1-1/r^2))./[1:l].^2 + (1-cond_nr/r^2)/(1-1/r^2);
    elseif strcmp(modeX0,'condition_control_linear')
        fct = @(l) ((cond_nr*r-1)/(r-1))+((1-cond_nr)/(r-1)).*[1:l];
    elseif strcmp(modeX0,'condition_control_log')
        fct = @(l) cond_nr.*exp(-log(cond_nr).*([1:l]-1)./(r-1));
    elseif strcmp(modeX0,'condition_control_log_plateau')
        fct = @(l) cond_nr.*[(exp(-log(cond_nr).*(((3/2).*[1:l/3])-1)./(r-1))),...
            (exp(-log(cond_nr).*((r/2).*ones(1,2*r/3-r/3-1)-1)./(r-1))),...
            (exp(-log(cond_nr).*([r/2:3/2:l]-1)./(r-1)))];
    elseif strcmp(modeX0,'condition_control_log_plateau_end') 
        fct = @(l) cond_nr.*[(exp(-log(cond_nr).*(((3/2).*[1:floor(2*l/3)])-1)./(r-1))),...
            (exp(-log(cond_nr).*(floor(r.*ones(1,r/3))-1)./(r-1)))];
    elseif strcmp(modeX0,'condition_control_log_plateau_NEW')
        fct = @(l) [cond_nr.*ones(1,r-7-8) logspace(log10(cond_nr),0,8) 2.*ones(1,l-8-(r-7-8))];            
    end
    S = fct(r);
    U0 = U.*S;
elseif strcmp(modeX0,'bivariate')
    h = @(x,y) 1./(1+(x-y).^2);
    M = zeros(d1,d2);
    for i=1:d1
        for j=1:d2
            M(i,j)=h(10.*(i-1)/d1,10.*(j-1)/d2);
        end
    end
    [U,S,V]=svd(M);
    U0 = U*S;
end
V0 = V;

end


