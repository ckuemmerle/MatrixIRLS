function [model, outs] = R3MC_rone_updates(prob,params)
% This function implements an iterative, rank-updating version of the 
% algorithm R3MC, see the section about "Rank updating" of the paper
% "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
% by Bamdev Mishra and Rodolphe Sepulchre, 53rd IEEE Conference on 
% Decision and Control. IEEE, 2014.
% Modified code of Bamdev Mishra.
% =========================================================================
% Author: Christian Kuemmerle, Johns Hopkins University, kuemmerle@jhu.edu,
% 2020.
% =========================================================================
rank_steps = 1;
max_rank = prob.r;
start_rank = params.start_rank;
ranklist = [];

for r_c = start_rank:rank_steps:max_rank
    prob.r      = r_c;
    if r_c == start_rank
        prob = initialize_R3MC(prob,params.random_initialization);
        outs = struct;
    end
    [model_c, outs_c] = R3MC_adp(prob, params);
    ranklist = [ranklist,r_c];
    [u, sig, v] = svds(model_c.M, 1);
    [prob.U0_init, prob.R0_init, prob.V0_init]...
        = svd_rank_1_update(model_c.U, model_c.R,...
        model_c.V, -sig*u, v);
    outs_c.X = outs_c.X.';
    if r_c == start_rank
        outs = outs_c;
    else
        outs_c.time = outs.time(end) + outs_c.time;
        outs = cell2struct(cellfun(@vertcat,struct2cell(outs),...
            struct2cell(outs_c),'uni',0),fieldnames(outs),1);
    end
end
outs.ranklist = ranklist;
outs.N_vec = outs.N;
outs.N = sum(outs.N_vec);
outs.X = outs.X.';
model = model_c;
end

