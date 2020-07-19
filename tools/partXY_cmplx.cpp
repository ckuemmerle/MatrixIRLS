/* -------------------------------------------------------------------------- */
/* partXY_mex mexFunction */
/* -------------------------------------------------------------------------- */

#include <mex.h>
#include <complex.h>
#include <tgmath.h>
#include <stdio.h>
/*#include "blas.h"*/

/* compute a part of X*Y */
void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin  [ ]
)
{
    double *Xtr, *Xti, *Yr, *Yi, *Z, *I, *J, *vr, *vi, LL;
//     std::complex<double> *Xt, *Y, *Z ,*v, *I, *J, LL;
//     double complex *I, *J, LL;
    ptrdiff_t m, n, r, L, p, ir, jr, k;
    ptrdiff_t inc = 1;

    if (nargin != 5 || nargout > 1)
        mexErrMsgTxt ("partXY(U', V', row, col, length(data))") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    
    Xtr = mxGetPr( pargin [0] );     // r x m
    Xti = mxGetPi( pargin [0] );     // r x m
    Yr  = mxGetPr( pargin [1] );     // r x n
    Yi  = mxGetPi( pargin [1] );     // r x n
    I  = mxGetPr( pargin [2] );     // row position
    J  = mxGetPr( pargin [3] );     // col position
    LL = mxGetScalar( pargin [4] ); // num of obvs
    L = (ptrdiff_t) LL;
    m  = mxGetN( pargin [0] );
    n  = mxGetN( pargin [1] );
    r  = mxGetM( pargin [0] ); 
    if ( r != mxGetM( pargin [1] ))
        mexErrMsgTxt ("rows of Xt must be equal to rows of Y") ;
    if ( r > m || r > n )
        mexErrMsgTxt ("rank must be r <= min(m,n)") ;
    
    /* ---------------------------------------------------------------- */
    /* output */
    /* ---------------------------------------------------------------- */

    pargout [0] = mxCreateDoubleMatrix(1, L, mxCOMPLEX);
    vr = mxGetPr( pargout [0] );
    vi = mxGetPi( pargout [0] );
    // v = mxGetPr( pargin [5] );
    
    /* C array indices start from 0 */
    for (p = 0; p < L; p++) {
        ir = ( I[p] - 1 ) * r;
        jr = ( J[p] - 1 ) * r;
         vr[p] = 0.0;
         vi[p] = 0.0;
        for (k = 0; k < r; k++) {
            vr[p] += Xtr[ ir + k ] * Yr[ jr + k ];
            vr[p] -= Xti[ ir + k ] * Yi[ jr + k ];
            vi[p] += Xtr[ ir + k ] * Yi[ jr + k ];
            vi[p] += Xti[ ir + k ] * Yr[ jr + k ];
        }
//             v[p] += Xt[ ir + k ] * Y[ jr + k ];
        /*v[p] = ddot(&r, Xt+ir, &inc, Y+jr, &inc);*/
    }
    
    return;
}

