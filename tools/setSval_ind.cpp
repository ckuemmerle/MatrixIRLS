/* -------------------------------------------------------------------------- */
/* setSval_ind mexFunction */
/* -------------------------------------------------------------------------- */

#include "mex.h"
#include "blas.h" 
#include "matrix.h"

/* set Write the values of v onto the indices of S given by ind */
void mexFunction(int nargout, mxArray *pargout [], int nargin, const mxArray *pargin [])
{
    double *Svalr;
    double *vr;
    double *ind;
    double  LL;  
    unsigned long L;
    unsigned long k;
    ptrdiff_t i;
    ptrdiff_t inc = 1;
    
    if (nargin != 4)
        mexErrMsgTxt ("Usage:  setSval_ind( S,ind, v, L ), set values of S(ind) to v") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    Svalr = mxGetPr( pargin [0] );
    ind   = mxGetPr( pargin [1] );
    vr    = mxGetPr( pargin [2] );
    
    LL    = mxGetScalar( pargin [3] );
    L = (unsigned long) LL;
    if (mxIsComplex(pargin[1])) {
        double *vi;
        double *Svali;
        vi    = mxGetPi( pargin [2] );
        if (!mxIsComplex(pargin[0])) {
            mexErrMsgIdAndTxt( "MATLAB:convec:inputsNotComplex",
              "For complex arithmetic, please also allocate complex memory in the sparse matrix.\n");
        }
        else {
            Svali = mxGetPi( pargin [0] );
            for (k = 0; k < L; k++) {
                i=ind[k]-1;
                Svalr[i] = vr[k]; 
                Svali[i] = vi[k];
            }
        }
    }
//         }
//         else {
//             Svali = mxCalloc(L,sizeof(double));
//         }
    else {
        for (k = 0; k < L; k++) {
            i=ind[k]-1;
            Svalr[i] = vr[k]; 
        }
    }
//     if (mxIsComplex(pargin[0])) {
//         Svali = mxGetPi( pargin [0] );
//     }
    
      
    
//     /* ---------------------------------------------------------------- */
//     /* output */
//     /* ---------------------------------------------------------------- */
//     if (mxIsComplex(pargin[1])) {
//         for (k = 0; k < L; k++) {
//             Svalr[k] = vr[k]; 
//             Svali[k] = vi[k];
//         }
//     }
//     else {
//         for (k = 0; k < L; k++) {
//             Svalr[k] = vr[k]; 
//         }
//     }
    
    return;
}
