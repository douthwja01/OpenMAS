/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_OMAS_updateGlobalStates_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

#ifndef _CODER_OMAS_UPDATEGLOBALSTATES_API_H
#define _CODER_OMAS_UPDATEGLOBALSTATES_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_OMAS_updateGlobalStates_api.h"

/* Type Definitions */
#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T

struct emxArray_boolean_T
{
  boolean_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_boolean_T*/

#ifndef typedef_emxArray_boolean_T
#define typedef_emxArray_boolean_T

typedef struct emxArray_boolean_T emxArray_boolean_T;

#endif                                 /*typedef_emxArray_boolean_T*/

#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T

struct emxArray_char_T
{
  char_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_char_T*/

#ifndef typedef_emxArray_char_T
#define typedef_emxArray_char_T

typedef struct emxArray_char_T emxArray_char_T;

#endif                                 /*typedef_emxArray_char_T*/

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

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  uint8_T objectID;
  emxArray_char_T *name;
  emxArray_char_T *class;
  uint8_T type;
  real32_T colour[3];
  emxArray_char_T *symbol;
  real_T radius;
  real_T detectionRadius;
  boolean_T idleStatus;
  real_T globalState[10];
  real_T R[9];
  emxArray_real_T *relativePositions;
  emxArray_boolean_T *objectStatus;
} struct1_T;

#endif                                 /*typedef_struct1_T*/

#ifndef typedef_emxArray_struct1_T
#define typedef_emxArray_struct1_T

typedef struct {
  struct1_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
} emxArray_struct1_T;

#endif                                 /*typedef_emxArray_struct1_T*/

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

#ifndef typedef_struct2_T
#define typedef_struct2_T

typedef struct {
  real_T duration;
  real_T dt;
  real_T startTime;
  real_T idleTimeOut;
  real_T endTime;
  emxArray_real_T *timeVector;
  real_T currentTime;
  real_T frequency;
  real_T numSteps;
  real_T currentStep;
} struct2_T;

#endif                                 /*typedef_struct2_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  emxArray_struct1_T *OBJECTS;
  struct2_T TIME;
  real_T conditionTolerance;
  emxArray_uint8_T *globalIDvector;
  boolean_T gui;
  boolean_T monteCarloMode;
  emxArray_char_T *outputPath;
  emxArray_char_T *phase;
  emxArray_char_T *systemFile;
  boolean_T threadPool;
  uint8_T totalAgents;
  uint8_T totalMiscs;
  uint8_T totalObjects;
  uint8_T totalObstacles;
  uint8_T totalWaypoints;
  uint8_T verbosity;
  real_T visabilityModifier;
  real_T warningDistance;
  emxArray_char_T *outputFile;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void OMAS_updateGlobalStates(struct0_T *SIM, uint8_T objectID, real_T
  globalVelocity_k[3], real_T quaternion_k[4], boolean_T idleStatus_k, struct1_T
  *METAObjUpdate);
extern void OMAS_updateGlobalStates_api(const mxArray * const prhs[5], int32_T
  nlhs, const mxArray *plhs[1]);
extern void OMAS_updateGlobalStates_atexit(void);
extern void OMAS_updateGlobalStates_initialize(void);
extern void OMAS_updateGlobalStates_terminate(void);
extern void OMAS_updateGlobalStates_xil_terminate(void);

#endif

/*
 * File trailer for _coder_OMAS_updateGlobalStates_api.h
 *
 * [EOF]
 */
