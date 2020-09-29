function prob = initialize_R3MC(prob,random_initialization)
% Provides the initialization for the algorithm decribed in the paper
% "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
% by Bamdev Mishra and Rodolphe Sepulchre, 53rd IEEE Conference on 
% Decision and Control. IEEE, 2014.
%
% =========================================================================
% Author: Christian Kuemmerle, Johns Hopkins University, kuemmerle@jhu.edu,
% 2020.
% =========================================================================
r = prob.r;
d1 = prob.d1;
d2 = prob.d2;

if ~random_initialization
    sparse_structure = sparse(prob.data_ls.rows, prob.data_ls.cols,...
        prob.data_ls.entries, d1, d2, prob.data_ls.nentries);
    [U,B,V] = svds(sparse_structure, r);
%     G = U*(B.^(0.5));
%     H = V*(B.^(0.5));
    prob.U0_init = U;
    prob.V0_init = V;
    prob.R0_init = B;
    fprintf('**** Initialization by taking %i dominant SVD\n', r);
else
    G = randn(d1, r); H = randn(d2, r);
    [Qg,Rg] = qr(G, 0);
    [Qh,Rh] = qr(H, 0);
    [q1,b,q2] = svd(Rg * Rh');
    prob.U0_init = Qg * q1;
    prob.V0_init = Qh * q2;
    prob.R0_init = b;
    fprintf('**** Random initialization\n');
end
end

