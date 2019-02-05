#include "mex.h"

#define RC_UP   0x0000004000
#define RC_NEAR 0x0000000000
#define RC_DOWN 0x0000002000

#define setRC(x) __asm__ volatile ("subl    $4, %%esp;       \
																	 movl     %%esp, %%ebx;    \
																	 stmxcsr  (%%ebx);         \
																	 andl     $0xffff9fff, (%%ebx); \
																	 orl      %[RC], (%%ebx);      \
																	 ldmxcsr  (%%ebx);             \
																	 addl     $4, %%esp;" \
																	 : : [RC] "n" (x))


void mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] ) {

  int rnd;

  rnd = mxGetScalar(prhs[0]);

  switch (rnd) {
    case  1 :  setRC(RC_UP);   break;
    case  0 :  setRC(RC_NEAR); break;
    case -1 :  setRC(RC_DOWN); break;
    default :  setRC(RC_NEAR); break;
  }

} /* setround */

