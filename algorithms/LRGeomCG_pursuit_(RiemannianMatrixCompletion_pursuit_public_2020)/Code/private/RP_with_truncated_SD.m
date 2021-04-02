function [xc,ithist] = RP_with_truncated_SD(prob, opts)
% Rank r approximation of steepest descent iteration
%
%   r = opts.rank_increase

t0 = tic();
rank_inc = opts.rank_increase;
ithist = prep_stats();

A_S = []; ind_i_S = []; ind_j_S = [];

% first iteration: xc = 0
xc.U = []; xc.S = []; xc.V = []; cur_rank = 0;
xc = prepx(prob, xc);

ithist = compute_stats(prob, t0, xc, ithist);    

opts_cg = opts;
opts_cg.verbosity = 0;

%%% WARNING: the xc.U and xc.V are not always orthonormal
itc = 1;
while cur_rank <= opts.max_rank
    d = xc.on_omega - prob.data;
    
        

    if opts.verbosity > 0
        fprintf(' Outer iteration: %3i. ', itc)    
        fprintf('Relative error = %d, \n', ithist.training_err(itc) )
    end
    
    itc = itc+1;
    if cur_rank + rank_inc > opts.max_rank
        break
    end
    
    % increase rank and make a search direction
    updateSval(prob.temp_omega, d, length(d));
    grad_fullspace = prob.temp_omega*1;

    [U,S,V] = compute_best_rank_approx(grad_fullspace,rank_inc,opts);
    if strcmpi(opts.search_strategy, '1D')
        tmin = exact_search_onlyTxM_LR(prob,xc,U*S,V);        

        xc.S = blkdiag(xc.S, tmin*S);
        xc.U = [xc.U U];
        xc.V = [xc.V V];        
    
    elseif strcmpi(opts.search_strategy, 'full_space')
        xc.U = [xc.U U];
        xc.V = [xc.V V];        
        [xc.S,A_S,ind_i_S,ind_j_S] = getoptS_full_space(prob,xc.U,xc.V,A_S,ind_i_S,ind_j_S);
        
    elseif strcmpi(opts.search_strategy, 'diagonal_space')
        xc.U = [xc.U U];
        xc.V = [xc.V V];       
        [xc.S,A_S] = getoptS_diagonal_space(prob,xc.U,xc.V,A_S);
    end

    cur_rank = cur_rank + rank_inc;
    xc = prepx(prob, xc);
    
    
    [xc,hist_cg,fail_cg,msg_cg] = LRGeomCG(prob,opts_cg,xc);
    
    
    
    ithist = compute_stats(prob, t0, xc, ithist);
        
end



end

function [U,S,V] = compute_best_rank_approx(M_omega,cur_rank,opts)
if opts.rand_pca
    [U,S,V] = pca(M_omega,cur_rank,true); 
else
    opts_svds.tol = 1e-10;
    [U,S,V,flag] = svds(M_omega, cur_rank, 'L', opts_svds);
    if flag, warning('svds did not converge'); end
end

end

function [S,A,ind_i,ind_j] = getoptS_full_space(prob,L,R,A,ind_i,ind_j)
% Can only be used if L and R are the same expect for the last column!!!      
    r = size(L,2); 
    if isempty(A)  
        ind_i = []; ind_j = [];
        A = zeros(length(prob.data),r^2); 
        ii = 1;
        for j = 1:r
            for i = 1:r
                ind_i = [ind_i i]; ind_j = [ind_j j];
                A(:,ii) = partXY(L(:,i)', R(:,j)', prob.Omega_i, prob.Omega_j, prob.m)';  
                ii = ii+1;
            end
        end
        
    else
        for i=1:r
            j = r;
            A = [A partXY(L(:,i)', R(:,j)', prob.Omega_i, prob.Omega_j, prob.m)'];
            ind_i = [ind_i i]; ind_j = [ind_j j];
        end
        for j=1:r-1
            i = r;
            A = [A partXY(L(:,i)', R(:,j)', prob.Omega_i, prob.Omega_j, prob.m)'];
            ind_i = [ind_i i]; ind_j = [ind_j j];
        end
        
    end
    
    
    %S = A\C ;  % seems very slow for some reason  -> can be done with QR
    %update (not in Matlab, see fortran code of Hammarling or Kressner)
    Svec = (A'*A)\(A'*prob.data);    
    S = zeros(r,r);    
    for ii=1:length(ind_i)
        S(ind_i(ii),ind_j(ii)) = Svec(ii);
    end
    
    
end



function [S,A] = getoptS_diagonal_space(prob,L,R,A_prev)
% Can only be used if L and R are the same expect for the last column!!!      
    if isempty(A_prev)
        r = size(L,2);
        A = zeros(length(prob.data),r);
        for i = 1:r
            d = partXY(L(:,i)', R(:,i)', prob.Omega_i, prob.Omega_j, prob.m)'; 
            A(:,i) = d(:);           
        end
    else
        A = [A_prev partXY(L(:,end)', R(:,end)', prob.Omega_i, prob.Omega_j, prob.m)'];        
    end

    
    %S = A\C ;  % seems very slow for some reason
    S = (A'*A)\(A'*prob.data);
    S = diag(S) ;
end


function ithist = prep_stats()    
    ithist.training_err = [];
    ithist.testing_err = [];
    ithist.time = [];
end

function ithist = compute_stats(prob, t0, xc, ithist)
    itc = length(ithist.time)+1;
    d = xc.on_omega - prob.data;
    ithist.training_err(itc) = 0.5*norm(d)/prob.norm_M_Omega;
    if prob.has_testing
        if isempty(xc.S)
            on_gamma = 0;
        else
            on_gamma = partXY(xc.S'*xc.U',xc.V', prob.Gamma_i, prob.Gamma_j, prob.m);
        end
        ithist.testing_err(itc) = 0.5*norm( on_gamma - prob.test_data' ) / norm(prob.test_data);
    else
        ithist.testing_err(itc) = Inf;
    end
    ithist.time(itc) = toc(t0);
end
