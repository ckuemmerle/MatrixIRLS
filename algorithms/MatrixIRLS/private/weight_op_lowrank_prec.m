function [H,Uprod,Vprod,dH,Hminc] = weight_op_lowrank_prec(U,V,S,lambda,expo,...
    type_mean,rho,U0,V0,varargin)
%(U,V,S,eps,expo,...
%    type_mean,rho,U0,V0)
%weight_op_lowrank_prec This precalculates some auxiliary matrices that
%are used in the functions weight_op_lowrank.m, given a harmonic mean or 
%Wasserstein mean weight matrix such that (U,V) are its left and 
%right singular vector matrices, 
%a reference rank-r matrix with singular vector matrices (U0,V0)%, and S
%the values for its Hadamard representation, used for the application of
%the matrix (lambda.*W + rho.*Id)^(expo).
%  Input:     U = (d1xR) matrix with the R principal left singular vectors
%                 of the iterate X^n
%             V = (d2xR) matrix with the R principal right singular vectors
%                 of the iterate X^n
%             S = (1x(R+1)) vector such that S(i) are the values that 
%                 define the Hadamard matrix H of the HM operator
%        lambda = Value for lambda (or eps in case of
%                       type_mean='2ndorder_modulus')
%          expo = parameter that determines to which exponent
%                 (lambda.*W + rho.*Id)^(expo) of lambda.*W + rho.*Id
%                 the operator corresponds  
%     type_mean = Parameter determine if the underlying weight matrix W is
%                 s.t.
%                 type_mean == 'harmonic': W is harmonic mean of one-sided
%                 weight matrices as in (Kuemmerle, Sigl 2017)
%                 type_mean == 'Wasserstein': W^(-1) is Wasserstein mean of
%                 the inverses of one-sided weights
%           rho = ADMM parameter rho that corresponds to the regularization
%                 (lambda.*W+rho.*Id).
%            U0 = (d1xR) matrix with the R principal left singular vectors
%                 of reference matrix X
%            V0 = (d2xR) matrix with the R principal right singular vectors
%                 of reference matrix X
% Output:     H = (1x4) cell such that
%                       H{1} corresponds to H_{R,R}
%                       H{2} corresponds to H_{R,R}-H_3
%                       H{3} corresponds to H_{R,R}-H_2
%                       H{4} corresponds to H_{R,R}-H_2-H_3+c*Id
%         Uprod = (RxR) matrix such that Uprod=U'*U0
%         Vprod = (RxR) matrix such that Vprod=V0'*V;
%            dH = (1x2) cell such that
%                       dH{1}(i) = H(R+1,i) for i=1..R+1,
%                       dH{2}(i) = H(R+1,i) - H(R+1,R+1) for i=1..R.
H = cell(1,4);
dH= cell(1,2);
R1 = size(U,2);
R2 = size(V,2);
Slen=length(S);

for k=1:4
    H{k}=zeros(R1,R2);
end
H2_tmp=zeros(R1,R2);
H3_tmp=zeros(R1,R2);

if nargin >= 8 
    if not(isempty(U0)) && not(isempty(V0))
        Uprod = U'*U0;
        Vprod = V0'*V;
    else
        Uprod = [];
        Vprod = [];
    end
else
    Uprod = [];
    Vprod = [];
end
if nargin >= 10
    qmean_para = varargin{1};
end
if nargin >= 11 && nargin < 14
    error('Not enough parameters for the anti-symmetric version.')
end
if nargin >= 14
    %%% In this case, the weights on the antisymmetric parts are increase
    %%% to match H_2^{(k)} as in the Ph.D. thesis [Kuemmerle 2019].
    increase_antisymmetricweights = 1;
    p=varargin{2};
    sing=varargin{3};
    eps=varargin{4};
    mode=varargin{5};
    if rho > 0
        error('Rho > 0 not supported for this type of weight matrices in the current implementation.')
    end        
else
    increase_antisymmetricweights = 0;
end
if nargin == 6
   rho = 0; 
end
if R1 > 0
    switch type_mean
    case 'harmonic'
        for i=1:R1
            for j=1:R2
                H{1}(i,j)   =(2*lambda/(S(i)+S(j))+rho)^(expo);
                H2_tmp(i,j) =(2*lambda/(S(j)+S(min(Slen,R1+1)))+rho)^(expo);
                H3_tmp(i,j) =(2*lambda/(S(i)+S(min(Slen,R2+1)))+rho)^(expo);
            end
        end
    case 'qmean'
        for i=1:R1
            for j=1:R2
                H{1}(i,j)   =(lambda*(2/(S(i)^(-qmean_para)+S(j)^(-qmean_para)))^(1/(-qmean_para))+rho)^(expo);
                H2_tmp(i,j) =(lambda*(2/(S(j)^(-qmean_para)+S(min(Slen,R1+1))^(-qmean_para)))^(1/(-qmean_para))+rho)^(expo);
                H3_tmp(i,j) =(lambda*(2/(S(i)^(-qmean_para)+S(min(Slen,R2+1))^(-qmean_para)))^(1/(-qmean_para))+rho)^(expo);
            end
        end
    case 'geometric'
        for i=1:R1
            for j=1:R2
                H{1}(i,j)   =(lambda/sqrt(S(i)*S(j))+rho)^(expo);
                H2_tmp(i,j) =(lambda/sqrt(S(j)*S(min(Slen,R1+1)))+rho)^(expo);
                H3_tmp(i,j) =(lambda/sqrt(S(i)*S(min(Slen,R2+1)))+rho)^(expo);
            end
        end
    case 'arithmetic'
        for i=1:R1
            for j=1:R2
                H{1}(i,j)   =(lambda*(1/S(i)+1/S(j))/2+rho)^(expo);
                H2_tmp(i,j) =(lambda*(1/S(min(Slen,R1+1))+1/S(j))/2+rho)^(expo);
                H3_tmp(i,j) =(lambda*(1/S(i)+1/S(min(Slen,R2+1)))/2+rho)^(expo);
            end
        end
    case 'Wasserstein'
        for i=1:R1
            for j=1:R2
                H{1}(i,j)   =(4*lambda/(S(i)+S(j)+2*sqrt(S(i)*S(j)))+rho)^(expo);
                H2_tmp(i,j) =(4*lambda/(S(j)+S(min(Slen,R1+1))+2*sqrt(S(min(Slen,R1+1))*S(j)))+rho)^(expo);
                H3_tmp(i,j) =(4*lambda/(S(i)+S(min(Slen,R2+1))+2*sqrt(S(i)*S(min(Slen,R2+1))))+rho)^(expo);
            end
        end
    case 'min'
        for i=1:R1
            for j=1:R2
                H{1}(i,j)   =(lambda/max(S(i),S(j))+rho)^(expo);
                H2_tmp(i,j) =(lambda/max(S(j),S(min(Slen,R1+1)))+rho)^(expo);
                H3_tmp(i,j) =(lambda/max(S(i),S(min(Slen,R2+1)))+rho)^(expo);
            end
        end
    case 'Hessianupper_1'
        for i=1:R1
            for j=1:R2
                H{1}(i,j)   =(lambda^(2-p)*(S(i)^(p-1)+S(j)^(p-1))/(S(i)+S(j)))^(expo);
                H2_tmp(i,j) =(lambda^(2-p)*(lambda^(p-1)+S(j)^(p-1))/(lambda+S(j)))^(expo);
                H3_tmp(i,j) =(lambda^(2-p)*(S(i)^(p-1)+lambda^(p-1))/(S(i)+lambda))^(expo);
            end
        end
    if increase_antisymmetricweights
        for i=1:R1
            for j=1:R2
                switch mode
                    case 'objective_thesis'
                        H2_tmp(i,j)    = (lambda*(sing(j)^(p-1)+eps^(p-1))/(sing(j)+eps))^(expo);
                    case 'pluseps'
                        H2_tmp(i,j)    = (lambda*((sing(j)+eps)^(p-1)+(2*eps)^(p-1))/(sing(j)+eps))^(expo);
                    case 'pluseps_squared'
                        H2_tmp(i,j)    = (lambda*(2*sing(j)/(sing(j)^2+eps^2)^(1-p/2)+2^(p/2)*eps^(p-1))/(sing(j)+eps))^(expo);
                    case 'JMLRpaper'
                        H2_tmp(i,j)    = (lambda*(2*sing(j)/(sing(j)^2+eps^2)^(1-p/2)+2^(p/2)*eps^(p-1))/(sing(j)+eps))^(expo);
                end
                if i > j
                    switch mode
                        case 'objective_thesis'
                            H{1}(i,j)    = (lambda*(sing(i)^(p-1)+sing(j)^(p-1))/(sing(i)+sing(j)))^(expo);
                        case 'pluseps'
                            H{1}(i,j)    = (lambda*((sing(i)+eps)^(p-1)+(sing(j)+eps)^(p-1))/(sing(i)+sing(j)))^(expo);
                        case 'pluseps_squared'
                            H{1}(i,j)    = (lambda*(2*sing(i)/(sing(i)^2+eps^2)^(1-p/2)+2*sing(j)/(sing(j)^2+eps^2)^(1-p/2))/(sing(i)+sing(j)))^(expo);
                        case 'JMLRpaper'
                            H{1}(i,j)    = (lambda*(2*sing(i)/(sing(i)^2+eps^2)^(1-p/2)+2*sing(j)/(sing(j)^2+eps^2)^(1-p/2))/(sing(i)+sing(j)))^(expo);
                    end
                end
            end
        end
    end
end
if R1 > 0 && R1==R2 
    R=R1;
    dH{1}=zeros(R+1,1);
    dH{2}=zeros(R,1);
    if  increase_antisymmetricweights
        dH{1}(1:R)=H3_tmp(1:R,1);
        dH{2}(1:R)=H2_tmp(1,1:R);
    else
        dH{1}(1:R)=H3_tmp(1:R,1);
        dH{1}(min(Slen,R+1))=(lambda/S(min(Slen,R+1))+rho)^(expo);
        dH{2}(1:R)=dH{1}(1:R)-dH{1}(min(Slen,R+1));
        H{2}=H{1}-H3_tmp;
        H{3}=H{1}-H2_tmp;
        H{4}=H{2}-H2_tmp+(lambda/S(min(Slen,R+1))+rho)^(expo).*ones(R,R);
        Hminc=H{1}-(lambda/S(min(Slen,R+1))+rho)^(expo).*ones(R,R);
    end
    if strcmp(type_mean,'2ndorder_modulus')
        cst_eps= (1/(lambda)^2)^(expo); 
        dH{1}(min(Slen,R+1))= cst_eps;
        dH{2}(1:R)=dH{1}(1:R)-dH{1}(min(Slen,R+1));
        H{4}=H{2}-H2_tmp+cst_eps.*ones(R,R);
        Hminc=H{1}-cst_eps.*ones(R,R);
    end
elseif R1 == 0
    dH{1} = (lambda/S(1)+rho)^(expo);
else
    error("Check dimensionalities of U and V, something is wrong.")
end
end

