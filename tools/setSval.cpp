/* -------------------------------------------------------------------------- */
/* SetSpv mexFunction */
/* -------------------------------------------------------------------------- */

#include "mex.h"
#include "blas.h" 
#include "matrix.h"

/* set values to sparse matrix S */
void mexFunction(int nargout, mxArray *pargout [], int nargin, const mxArray *pargin [])
{
    double *Svalr;
    double *vr;
    double  LL;  
    unsigned long L;
    unsigned long k;
    
    if (nargin != 3)
        mexErrMsgTxt ("Usage:  SetSpv ( S, v, L ) from v to S") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    Svalr = mxGetPr( pargin [0] );
    vr    = mxGetPr( pargin [1] );
    LL    = mxGetScalar( pargin [2] );
    L = (unsigned long) LL;
    if (mxIsComplex(pargin[1])) {
        double *vi;
        double *Svali;
        vi    = mxGetPi( pargin [1] );
        if (!mxIsComplex(pargin[0])) {
            mexErrMsgIdAndTxt( "MATLAB:convec:inputsNotComplex",
              "For complex arithmetic, please also allocate complex memory in the sparse matrix.\n");
        }
        else {
            Svali = mxGetPi( pargin [0] );
            for (k = 0; k < L; k++) {
                Svalr[k] = vr[k]; 
                Svali[k] = vi[k];
            }
        }
    }
//         }
//         else {
//             Svali = mxCalloc(L,sizeof(double));
//         }
    else {
        for (k = 0; k < L; k++) {
            Svalr[k] = vr[k]; 
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
