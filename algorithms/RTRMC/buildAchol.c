/*=================================================================
% function Achol = buildAchol(Chat, U, sqlambda, I, J, n)
%
% Specific function to reduce the computation time in a bottleneck
% portion of the Matlab code for our low-rank matrix completion
% software RTRMC.
%
% Chat is a real, double >= 0 col vector, U is an orthonormal (real,double)
% m-by-r matrix, sqlambda is a positive real number, I and J are uint32
% vectors of the same size as Chat; together, Chat, I and J define an
% m-by-n sparse matrix; n must be a uint32.
%
% I and J are assumed to be sorted just like the output of the FIND
% function in Matlab.
% 
% Achol is a cell containing n upper triangular r-by-r matrices: the
% Cholesky factors of the diagonal blocks of A. See the notes
% accompanying the software for details.
%
% Compile with: mex -lmwlapack -lmwblas -largeArrayDims buildAchol.c
%
% May 19, 2011  Nicolas Boumal, UCLouvain
 *=================================================================*/

#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "lapack.h"
#include "blas.h"

/* Input Arguments */

#define	pChat	    prhs[0]
#define	pU          prhs[1]
#define	psqlambda   prhs[2]
#define pI          prhs[3]
#define pJ          prhs[4]
#define pn          prhs[5]

/* Output Arguments */

#define pAchol  	plhs[0]

/* Helper function */
/* mwSize max(mwSize a, mwSize b) { return a > b ? a : b; } */
#define max(a, b) ((a)>(b)?(a):(b))


void mexFunction(
          int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray* prhs[] )
{
    /* Counters and sizes*/
    mwSize m, n, r, known, start, end;
    mwIndex i, j, k;
    
    /* Data arrays */
    double *Chat, *U;
    uint32_T *I, *J, *nn;
    double sqlambda;
    
    /* Temporary computation variables */
    mxArray *pscaledUi, *pAi;
    double *scaledUi, *Ai;
    mwSize nb;
    double coeff;
    
    /* BLAS/LAPACK stuff*/
    ptrdiff_t N, M, R, NB, info = 0;
    double one = 1.0;
    double zero = 0.0;
    char *uplo = "U";
    char *TRANS = "T";
    char *NOTRANS = "N";
    
    
    if(nrhs != 6 || nlhs != 1)
        mexErrMsgTxt("Invalid number of input/output parameter. Need 6 inputs and 1 output.");
    
    
    r = mxGetN(pU);
    m = mxGetM(pU);
    nn = (uint32_T*) mxGetData(pn);
    n = (mwSize) (nn[0]);
    
    known = max(mxGetM(pI), mxGetN(pI));
    if(max(mxGetM(pJ), mxGetN(pJ)) != known || max(mxGetM(pChat), mxGetN(pChat)) != known)
        mexErrMsgTxt("I, J and Chat must be vectors of the same length.");
    
    
    sqlambda = mxGetScalar(psqlambda);
    
    Chat = mxGetPr(pChat);
    U = mxGetPr(pU);
    I = (uint32_T*) mxGetData(pI);
    J = (uint32_T*) mxGetData(pJ);
    
    
    N = (ptrdiff_t) n;
    M = (ptrdiff_t) m;
    R = (ptrdiff_t) r;
    
    
    /* dummy matrix, for memory allocation */
    pscaledUi = mxCreateDoubleMatrix(m, r, mxREAL);
    scaledUi  = mxGetPr(pscaledUi);
    
    /* Create a suitable cell */
    pAchol = mxCreateCellMatrix(n, 1);
    for(i = 0; i < n; ++i)
        mxSetCell(pAchol, i, NULL);
    
    start = 0;
    end = 0;
    for(i = 0; i < n; ++i)
    {
        /* allocate space for an r-by-r matrix in each cell entry */
        mxSetCell(pAchol, i, mxCreateDoubleMatrix(r, r, mxREAL));
        pAi = mxGetCell(pAchol, i);
        Ai  = mxGetPr(pAi);
        
        /* Compute the entries of the i-th diagonal block of A
           nb is the number of nonzero entries in the i-th column of Chat
           seen as a sparse matrix of size m-by-n with nonzero entries I, J. */
        
        /* TODO TODO TODO
           It would be much better to have an additionnal guard element at
           the end of the J vector ... :(
           Actually, it would be even better to compute the Jc vector
           like the one you get in the Matlab sparse representation.
           And even better would be to precompute it just once. */
        
        while(end < known && J[end] == i+1) {
            ++end;
        }
        
        nb = end-start;
        
        if(nb > 0)
        {
/*///////////////////////////////////////////////////////////////////////*/
            for(j = 0; j < nb; ++j)
            {
                coeff = sqrt(Chat[start+j]);
                for(k = 0; k < r; ++k)
                {
                    scaledUi[j+k*m] = coeff * U[I[start+j]-1+k*m];
                }
            }
/*///////////////////////////////////////////////////////////////////////*/
            
            /* compute the matrix product of the interesting part of
               scaledUi with itself (transposed): Ai = scaledUi.'*scaledUi */
            NB = (ptrdiff_t) nb;
            dgemm(TRANS, NOTRANS, &R, &R, &NB, &one,
                  scaledUi, &M, scaledUi, &M, &zero, Ai, &R);
        }
        else
        {
            for(j = 0; j < r*r; ++j)
                Ai[j] = 0.0;
        }
        
        start = end;
        
        /* add sqlambda*eye(r) to Ai */
        for(j = 0; j < r; ++j)
            Ai[j*(r+1)] += sqlambda;
        
        /* remplace Ai by its Cholesky factor by calling Lapack's dpotrf. */
        dpotrf(uplo, &R, Ai, &R, &info);
        
        /* write zeroes on the lower triangular parts of Ai (optional) */
        for(j = 1; j < r; ++j)
            for(k = 0; k < j; ++k)
                Ai[j+r*k] = 0.0;
        
        if(info < 0)
            mexErrMsgTxt("dpotrf (in buildAchol): invalid argument.");
        if(info > 0)
            mexErrMsgTxt("dpotrf (in buildAchol): a leading minor is not positive definite.");
    }
    
    mxDestroyArray(pscaledUi);
    
    return;
}


