/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_OMAS_sphereTriangleIntersection_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 07-Nov-2018 17:03:23
 */

#ifndef _CODER_OMAS_SPHERETRIANGLEINTERSECTION_API_H
#define _CODER_OMAS_SPHERETRIANGLEINTERSECTION_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_OMAS_sphereTriangleIntersection_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern boolean_T OMAS_sphereTriangleIntersection(real_T P[3], real_T r, real_T
  A[3], real_T B[3], real_T C[3]);
extern void OMAS_sphereTriangleIntersection_api(const mxArray * const prhs[5],
  int32_T nlhs, const mxArray *plhs[1]);
extern void OMAS_sphereTriangleIntersection_atexit(void);
extern void OMAS_sphereTriangleIntersection_initialize(void);
extern void OMAS_sphereTriangleIntersection_terminate(void);
extern void OMAS_sphereTriangleIntersection_xil_terminate(void);

#endif

/*
 * File trailer for _coder_OMAS_sphereTriangleIntersection_api.h
 *
 * [EOF]
 */
