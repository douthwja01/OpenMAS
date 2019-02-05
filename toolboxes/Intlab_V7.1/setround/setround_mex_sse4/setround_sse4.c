#include "mex.h"
#include <float.h>
#include <omp.h>
#pragma STDC FENV_ACCESS ON


void mexFunction(int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[])
{

  int rnd, err, n, i;
  unsigned int control_word;

  rnd = mxGetScalar(prhs[0]);
  n=omp_get_max_threads();
  //printf("%d\n",n);
  
  #pragma omp parallel
  {
  #pragma omp for private(i)
  for (i=0;i<n;i++)
  {
      //printf("OK\n");
      switch (rnd) {
          case -1 :
              err = _controlfp_s(&control_word, _RC_DOWN, _MCW_RC);
              break;
          case  0 :
              err = _controlfp_s(&control_word, _RC_NEAR, _MCW_RC);
              break;
          case  1 :
              err = _controlfp_s(&control_word, _RC_UP, _MCW_RC);
              break;
          case  2 :
              err = _controlfp_s(&control_word, _RC_CHOP, _MCW_RC);
              break;
          default :
              err = _controlfp_s(&control_word, _RC_NEAR, _MCW_RC);
              break;
      }
  }
  }

} /* setround */

