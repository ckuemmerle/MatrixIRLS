function err = comp_error(M,X,varargin)

%computes relative frobenius norm between M (exact) and X using Lanczos SVD
%or traces (default)

if isfield(X,'U') && isfield(X,'S') && isfield(X,'V') %OptSpace,SET,SVT
    X.left = X.U*X.S;
    X.right = X.V;
    clear X.U X.S X.V  
elseif isfield(X,'Q') && isfield(X,'R') %GROUSE
    X.left = X.Q;
    X.right = X.R; 
    clear X.Q X.R
end

if isfield(X,'left') && isfield(X,'right') %OptSpace,SET,SVT,EALM,IALM,GROUSE
    if isfield(M,'left') && isfield(M,'right')
        if nargin == 3 && strcmpi(varargin{1},'exact') 
            err = exact_err(X,M);
        else
            err = fast_err(X,M);
        end
        
    elseif isfield(M,'exact')
        err = norm(M.exact-X.left*X.right','fro')/norm(M.exact,'fro'); 
        
    elseif isfield(M, 'fun')
        % we only have a function for M, sample the error
        n_1 = size(X.left,1); n_2 = size(X.left,1);
        m = min(100*n_1,n_1*n_2);
        err_ind_Omega = randsampling(n_1,n_2,m)';
        err_ind_Omega = err_ind_Omega(1:m); 
        [row, col] = ind2sub([n_1,n_2], err_ind_Omega);
        
        X_omega = partXY(X.left', X.right', row, col, m); 
        xx = linspace(0,1,n_1);
        M_omega = M.fun(xx(row),xx(col));
        err = norm( X_omega - M_omega )/norm( M_omega );
        

    end
elseif isfield(X,'full') %cvx, APGL
    if isfield(M,'exact')
        err = norm(M.exact-X.full,'fro')/norm(M.exact,'fro');
     elseif isfield(M,'left') && isfield(M,'right')
         M.exact = M.left*M.right';
         err = norm(M.exact-X.full,'fro')/norm(M.exact,'fro'); 
    end
end

end


function err = exact_err(X,M)
[Ul,Sl,Vl] = svd([X.left M.left],0);
[Ur,Sr,Vr] = svd([X.right -M.right],0);
nrm = norm(Sl*Vl'*Vr*Sr,'fro');
normMsqr = sum(sum((M.left'*M.left).*(M.right'*M.right)));
err = nrm/sqrt(normMsqr);
end

% function err = exact_err_lansvd(X,M) %based on partial SVD
% %due to lansvd
% clear global; 
% global N Y;
% N = M;
% Y = X; 
% opt.tol = 1e-12;
% opt.lanmax = round(1.2*(size(M.left,2)+size(X.left,2)))+10;
% [~,S,~,bnd] = lansvd('Bmap','BTmap',M.size1,M.size2,size(M.left,2)+size(X.left,2),'L',opt);
% norm(bnd);
% sigma = diag(S);
% nrm = sum(sigma.*sigma);
% normMsqr = sum(sum((M.left'*M.left).*(M.right'*M.right)));
% err = sqrt(nrm/normMsqr);
% end

function omega = randsampling(n_1,n_2,m)

omega = ceil(rand(m, 1) * n_1 * n_2);
omega = unique(omega);
while length(omega) < m    
    omega = [omega; ceil(rand(m-length(omega), 1)*n_1*n_2);];
    omega = unique(omega);
end
end
