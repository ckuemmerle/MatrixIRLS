/*=================================================================
% function X = cholsolve(R, B)
% Computes X such that A{i}X(:,i) = B(:,i), where R{i} is the upper
% triangular Cholesky factor of the real symmetric positive-definite
% matrix A{i}, i.e., A{i} = R{i}.'*R{i}.
% R is a cell array of m upper cholesky factors of real, symmetric,
% positive-definite matrices A{i}, i=1..m, of size n-by-n;
% B and X are n-by-m real matrices.
%
% Compile with: mex -lmwlapack -largeArrayDims cholsolvecell.c
%
% Feb. 17, 2011  Nicolas Boumal, UCLouvain
 *=================================================================*/

#include "mex.h"
#include "string.h"
#include "matrix.h"
#include "lapack.h"

/* Input Arguments */

#define	R	prhs[0]
#define	B	prhs[1]

/* Output Arguments */

#define	X	plhs[0]

void mexFunction(
          int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray* prhs[] )
{
    mwIndex k1, k2; /* debug */
    double *b;
    
    mwIndex i;
    mwSize m, n;
    mxArray *Ri;
    double *Rvals, *Bvals, *Xvals;
    char *uplo = "U";

    /* !! very important to use the proper typedef since confusion between
       32 and 64 bits for integers is quite common. */
    ptrdiff_t info = 0;
    ptrdiff_t one = 1;
	ptrdiff_t N = 0;
    
    /* Check for proper number of arguments */
    if (nrhs != 2) { 
        mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs != 1) {
        mexErrMsgTxt("Single output argument required."); 
    }
    
    if(!mxIsCell(R))
        mexErrMsgTxt("CHOLSOLVECELL: R must be a cell containing m upper triangular matrices of size n-be-n.");
        
    /* retrieve dimensions of the matrices */
    n = mxGetM(B);
    m = mxGetN(B);
    if(m <= 0 || n <= 0 || mxGetNumberOfElements(R) != m)
        mexErrMsgTxt("CHOLSOLVECELL: inconsistent matrix / cell dimensions.");
    
    /* dpotrs will store the result in the B matrix we pass to it;
       hence we copy B in a new matrix X and pass X to dpotrs.
       As a result, B will be unaltered and dpotrs will write the
       system solutions in the output matrix directly. */
    X = mxCreateDoubleMatrix(n, m, mxREAL);
    Xvals = mxGetPr(X);
    Bvals = mxGetPr(B);
    memcpy(Xvals, Bvals, m*n*sizeof(double));
    
    /* Call LAPACK subroutine (this is where the magic happens) */
	N = (ptrdiff_t) n;
    for(i = 0; i < m; ++i)
    {
        Ri = mxGetCell(R, i);
        if(mxGetM(Ri) != n || mxGetN(Ri) != n)
            mexErrMsgTxt("CHOLSOLVECELL: inconsistent matrix / cell dimensions.");
        Rvals = mxGetPr(Ri);
        b = Xvals+n*i;
        dpotrs(uplo, &N, &one, Rvals, &N, b, &N, &info);
        if(info < 0)
            mexErrMsgTxt("dpotrs (in cholsolvecell): invalid argument.");
    }
    
    return;
}


