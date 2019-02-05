#include <float.h>              /* Switch rounding mode, Dirk Husung */
#include "mex.h"

extern void setRC (unsigned short rc);
#pragma aux setRC =              \
  "sub   esp, 2"                 \
  "mov   ebx, esp"               \
  "fstcw word ptr [ebx]"         \
  "and   word ptr [ebx], 0f3ffh" \
  "or    word ptr [ebx], ax"     \
  "fldcw word ptr [ebx]"         \
  "pop   ax"                     \
  parm   [ax]                    \
  modify [ebx];

void mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] ) {

  int rnd;

  rnd = mxGetScalar(prhs[0]);

  switch (rnd) {
    case  1 :  setRC (RC_UP);   break;
    case  0 :  setRC (RC_NEAR); break;
    case -1 :  setRC (RC_DOWN); break;
    default :  setRC (RC_NEAR); break;
  }

} /* setround */

