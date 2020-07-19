/*=================================================================
% function R = spbuildmatrix(Ct, U, lambda)
%
% Specific function to reduce the computation time in a bottleneck
% portion of the Matlab code for our low-rank matrix completion
% software.
%
% Ct is a *sparse* m-by-n real matrix, U is an orthonormal (real)
% m-by-r matrix and lambda is a strictly positive real number.
% R is a cell containing n upper triangular r-by-r matrices: the
% Cholesky factors of the diagonal blocks of A.'*A. See the notes
% accompanying the software for details.
%
% The equivalent Matlab code is the following:
%
%         AtA = cell(n, 1);
%         for i = 1 : n
%             
%             I = find(Ct(:, i));
%             scaledUi = (Ct(I, i)*ones(1, r)) .* U(I, :);
%             AtAi = scaledUi.'*scaledUi + lambda*eye(r);
%             
%             % We actually store the Cholesky factorization
%             AtA{i} = chol(AtAi);
%             
%         end
%
% 
% Compile with: mex -lmwlapack -lmwblas -largeArrayDims spbuildmatrix.c
%
% Feb. 17, 2011  Nicolas Boumal, UCLouvain
% Bug correction on Oct. 15, 2014, for the case of non-uniform Ct.
 *=================================================================*/

#include "mex.h"
#include "matrix.h"
#include "lapack.h"
#include "blas.h"

/* Input Arguments */

#define	Ct	    prhs[0]
#define	U	    prhs[1]
#define	lambda	prhs[2]

/* Output Arguments */

#define R	plhs[0]

void mexFunction(
          int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray* prhs[] )
{
    mwIndex i, j, k;
    mwSize m, n, r;
    mxArray *Ri;
    double *Rivals, *Ctvals, *Uvals;
    char *uplo = "U";
    char *TRANS = "T";
    char *NOTRANS = "N";
    ptrdiff_t N, NB, M, RR, info = 0;
    double lambdaval;
    mwIndex *CtIr, *CtJc;
    mxArray *scaledUi;
    double *scaledUivals;
    double one = 1.0;
    double zero = 0.0;
    mwSize nb;
    mwIndex *Ctcol;
    double coeff;
    
    if(nrhs != 3 || nlhs != 1)
        mexErrMsgTxt("Invalid number of input/output parameter. Need 3 inputs and 1 output.");
    
    m = mxGetM(Ct);
    n = mxGetN(Ct);
    r = mxGetN(U);
    if(mxGetM(U) != m)
        mexErrMsgTxt("Dimensions mismatch in input arguments.");
    
    if(!mxIsSparse(Ct))
        mexErrMsgTxt("Ct ought to be a sparse matrix.");
    
    lambdaval = mxGetScalar(lambda);
    
    CtIr = mxGetIr(Ct);
    CtJc = mxGetJc(Ct);
    Ctvals = mxGetPr(Ct);
    Uvals = mxGetPr(U);
    
    N  = (ptrdiff_t) n;
    M  = (ptrdiff_t) m;
    RR = (ptrdiff_t) r;
    
    /* dummy matrix, for memory allocation */
    scaledUi = mxCreateDoubleMatrix(m, r, mxREAL);
    scaledUivals = mxGetPr(scaledUi);
    
    /* Create a suitable cell */
    R = mxCreateCellMatrix(n, 1);
    for(i = 0; i < n; ++i)
        mxSetCell(R, i, NULL);
    
    for(i = 0; i < n; ++i)
    {
        /* allocate space for an r-by-r matrix in each cell entry */
        mxSetCell(R, i, mxCreateDoubleMatrix(r, r, mxREAL));
        Ri = mxGetCell(R, i);
        Rivals = mxGetPr(Ri);
        
        /* compute the entries of the matrix R{i}
           nb is the number of nonzero entries in the i-th column of Ct
           Ctcol is a pointer to an array of nb elements: the indices of
           rows where a nonzero element exists on column i of Ct. */
        nb = CtJc[i+1]-CtJc[i];
        /* mexPrintf("%d ", nb); */
        Ctcol = CtIr+CtJc[i];
        if(nb > 0)
        {
            for(j = 0; j < nb; ++j)
            {
                /* coeff = Ct(Ctcol[j], i) */
                /* Correction on Oct. 15, 2014 */
                coeff = Ctvals[CtJc[i]+j];
                /* mexPrintf("%g ", coeff); */
                for(k = 0; k < r; ++k)
                {
                    /* scaledUi(j, k) = coeff * U(Ctcol[j], k) */
                    scaledUivals[j+k*m] = coeff * Uvals[Ctcol[j]+k*m];
                }
            }
            /*   mexPrintf("\n"); */
            /* compute the matrix product of the interesting part of
               scaledUi with itself (transposed): Ri = scaledUi.'*scaledUi */
            NB = (ptrdiff_t) nb;
            dgemm(TRANS, NOTRANS, &RR, &RR, &NB, &one,
                  scaledUivals, &M, scaledUivals, &M, &zero, Rivals, &RR);
        }
        else
        {
            for(j = 0; j < r*r; ++j)
                Rivals[j] = 0.0;
        }
        
        /* add lambda*eye(r) to Ri */
        for(j = 0; j < r; ++j)
            Rivals[j*(r+1)] += lambdaval;
        
        /* remplace R{i} by its Cholesky factor by calling Lapack's dpotrf. */
        dpotrf(uplo, &RR, Rivals, &RR, &info);
        /* write zeroes on the lower triangular parts of Ri */
        for(j = 1; j < r; ++j)
            for(k = 0; k < j; ++k)
                Rivals[j+r*k] = 0.0;
        
        if(info < 0)
            mexErrMsgTxt("dpotrf (in spbuildmatrix): invalid argument.");
        if(info > 0)
            mexErrMsgTxt("dpotrf (in spbuildmatrix): a leading minor is not positive definite.");
    }
    
    mxDestroyArray(scaledUi);
    
    return;
}


