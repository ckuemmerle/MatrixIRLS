function [Phi,Omega] = ...
    sample_phi_MatrixCompletion(d1,d2,m,mode,varargin)
%sample_phi_MatrixCompletion. This function randomly samples a (d1 x d2)
%sparse matrix with ones at m randomly chosen coordinates (uniform without
%replacement). If the 'resample' mode is chosen, this sampling procedure is
%repeated until Phi has at least r non-zero entries in each row and each
%column, where r is a specified positive integer, but this resampling
%procedure is repeated at most max_nr_resample times.

%           Input:
%           d1      = number of rows of Phi
%           d2      = number of columns of Phi 
%           m       = number of nonzero entries of Phi
%           mode    = 'resample', if resampling procedure is preferred, 
%                   = [], else.
%           If mode == 'resample':
%                   r = varargin{1}, lower bound on entries per column and
%                                    row
%     max_nr_resample = varargin{2}, upper bound on number of resamplings.
%
%           Output:
%           Phi     = (d1 x d2) sparse matrix: completion mask with 
%                       ones at "seen" entries and
%                       zeros at "unseen" entries.
%         Omega     = (m x 1) vector with linear indices of non-zero
%                       entries of Phi.

if strcmp(mode,'resample')
    r                = varargin{1};
    max_nr_resample  = varargin{2};
   
   
    rejec_counter=0;
    reject=1;
    while (reject==1)
        Omega = (sort(randperm(d1*d2,m)))';
        i_Omega=mod(Omega,d1);
        i_Omega(i_Omega==0)=d1;
        j_Omega=floor((Omega-1)/d1)+1;
        Phi = sparse(i_Omega,j_Omega,ones(m,1),d1,d2);
        nr_entr_col = sum(Phi,1)';
        nr_entr_row = sum(Phi,2);

        if (isempty(find(nr_entr_row<r+1,1)) == 0) || (isempty(find(nr_entr_col<r+1,1)) == 0)
            rejec_counter=rejec_counter+1;
        else
            reject=0;
        end
        if rejec_counter >= max_nr_resample
            break
        end
    end
    if rejec_counter >= 1
        disp(['Rejection counter of sampled matrix completion masks with',num2str(m),' entries: ',num2str(rejec_counter)]);
    end
else
    Omega = sort(randperm(d1*d2,m));
    i_Omega=mod(Omega,d1);
    i_Omega(i_Omega==0)=d1;
    j_Omega=floor((Omega-1)/d1)+1;
    Phi = sparse(i_Omega,j_Omega,ones(m,1),d1,d2);
end

end
