/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_OMAS_sphereTriangleIntersection_mex.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 07-Nov-2018 17:03:23
 */

/* Include Files */
#include "_coder_OMAS_sphereTriangleIntersection_api.h"
#include "_coder_OMAS_sphereTriangleIntersection_mex.h"

/* Function Declarations */
static void c_OMAS_sphereTriangleIntersecti(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[5]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[1]
 *                int32_T nrhs
 *                const mxArray *prhs[5]
 * Return Type  : void
 */
static void c_OMAS_sphereTriangleIntersecti(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[5])
{
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4,
                        31, "OMAS_sphereTriangleIntersection");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 31,
                        "OMAS_sphereTriangleIntersection");
  }

  /* Call the function. */
  OMAS_sphereTriangleIntersection_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  OMAS_sphereTriangleIntersection_terminate();
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(OMAS_sphereTriangleIntersection_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  OMAS_sphereTriangleIntersection_initialize();

  /* Dispatch the entry-point. */
  c_OMAS_sphereTriangleIntersecti(nlhs, plhs, nrhs, prhs);
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_OMAS_sphereTriangleIntersection_mex.c
 *
 * [EOF]
 */
