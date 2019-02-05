#include <omp.h>
#include <fenv.h>
#include "mex.h"

#pragma STDC FENV_ACCESS ON

void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    int rnd[] = {FE_DOWNWARD,FE_TONEAREST,FE_UPWARD,FE_TOWARDZERO}; 
    register int mode; 
   
    /*mode = *(mxGetPr(prhs[0]));*/
    mode = (int) mxGetScalar(prhs[0]);
    mode = mode + 1;

#pragma omp parallel shared(mode,rnd) 
{
    fesetround(rnd[mode]);
}

} 

/* the same like as fesetround, but in assembler

   __asm__ ("fnstcw %0" : "=m" (cw));
    cw &= 0x3ff;
    cw |= mode;
    __asm__ ("fldcw %0" : : "m" (cw));
    __asm__ ("stmxcsr %0" : "=m" (mxcsr));
    mxcsr &= ~ 0x6000;
    mxcsr |= mode << 3;
    __asm__ ("ldmxcsr %0" : : "m" (mxcsr));
 */

