/* -------------------------------------------------------------------------- */
/* SetSpv mexFunction */
/* -------------------------------------------------------------------------- */

#include "mex.h"
#include "blas.h" 
#include "matrix.h"

/* set the values of matrix X that correspond to indices in "ind"
 * the ones of Z. Index set of size L.  */
void mexFunction(int nargout, mxArray *pargout [], int nargin, const mxArray *pargin [])
{
    double *Xr;
    double *Zr;
    double *ind;
    double  LL;  
    unsigned long L;
    unsigned long k;
    ptrdiff_t i;
    ptrdiff_t inc = 1;
    
    if (nargin != 4)
        mexErrMsgTxt ("Usage:  setindXtoZ(X,Z,ind,L)") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    Xr   = mxGetPr( pargin [0] );
    Zr   = mxGetPr( pargin [1] );
    ind = mxGetPr( pargin [2] );
    LL  = mxGetScalar( pargin [3] );
    L = (unsigned long) LL;
    if (mxIsComplex(pargin[1])) {
        double *Xi;
        double *Zi;
        Zi    = mxGetPi( pargin [1] );
        if (!mxIsComplex(pargin[0])) {
            mexErrMsgIdAndTxt( "MATLAB:convec:inputsNotComplex",
              "For complex arithmetic, please also allocate complex memory in the sparse matrix.\n");
        }
        else {
            Xi = mxGetPi( pargin [0] );
            for (k = 0; k < L; k++) {
                i=ind[k]-1;
                Xr[i] = Zr[i]; 
                Xi[i] = Zi[i];
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
            Xr[i] = Zr[i]; 
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
