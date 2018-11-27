/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_OMAS_getTrajectoryData_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Nov-2018 17:59:35
 */

#ifndef _CODER_OMAS_GETTRAJECTORYDATA_API_H
#define _CODER_OMAS_GETTRAJECTORYDATA_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_OMAS_getTrajectoryData_api.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

#ifndef struct_emxArray_uint8_T
#define struct_emxArray_uint8_T

struct emxArray_uint8_T
{
  uint8_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_uint8_T*/

#ifndef typedef_emxArray_uint8_T
#define typedef_emxArray_uint8_T

typedef struct emxArray_uint8_T emxArray_uint8_T;

#endif                                 /*typedef_emxArray_uint8_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void OMAS_getTrajectoryData(emxArray_real_T *trajectoryMatrix,
  emxArray_uint8_T *globalIDvector, uint8_T objectID, real32_T step,
  emxArray_real_T *objectStateData, uint32_T *indexOfLastState);
extern void OMAS_getTrajectoryData_api(const mxArray * const prhs[4], int32_T
  nlhs, const mxArray *plhs[2]);
extern void OMAS_getTrajectoryData_atexit(void);
extern void OMAS_getTrajectoryData_initialize(void);
extern void OMAS_getTrajectoryData_terminate(void);
extern void OMAS_getTrajectoryData_xil_terminate(void);

#endif

/*
 * File trailer for _coder_OMAS_getTrajectoryData_api.h
 *
 * [EOF]
 */
