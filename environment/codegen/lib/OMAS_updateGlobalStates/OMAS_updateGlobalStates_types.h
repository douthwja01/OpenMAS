/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_updateGlobalStates_types.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

#ifndef OMAS_UPDATEGLOBALSTATES_TYPES_H
#define OMAS_UPDATEGLOBALSTATES_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T

struct emxArray_boolean_T
{
  boolean_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
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
  char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
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
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
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
  unsigned char objectID;
  emxArray_char_T *name;
  emxArray_char_T *class;
  unsigned char type;
  float colour[3];
  emxArray_char_T *symbol;
  double radius;
  double detectionRadius;
  boolean_T idleStatus;
  double globalState[10];
  double R[9];
  emxArray_real_T *relativePositions;
  emxArray_boolean_T *objectStatus;
} struct1_T;

#endif                                 /*typedef_struct1_T*/

#ifndef typedef_emxArray_struct1_T
#define typedef_emxArray_struct1_T

typedef struct {
  struct1_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
} emxArray_struct1_T;

#endif                                 /*typedef_emxArray_struct1_T*/

#ifndef struct_emxArray_uint8_T
#define struct_emxArray_uint8_T

struct emxArray_uint8_T
{
  unsigned char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
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
  double duration;
  double dt;
  double startTime;
  double idleTimeOut;
  double endTime;
  emxArray_real_T *timeVector;
  double currentTime;
  double frequency;
  double numSteps;
  double currentStep;
} struct2_T;

#endif                                 /*typedef_struct2_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  emxArray_struct1_T *OBJECTS;
  struct2_T TIME;
  double conditionTolerance;
  emxArray_uint8_T *globalIDvector;
  boolean_T gui;
  boolean_T monteCarloMode;
  emxArray_char_T *outputPath;
  emxArray_char_T *phase;
  emxArray_char_T *systemFile;
  boolean_T threadPool;
  unsigned char totalAgents;
  unsigned char totalMiscs;
  unsigned char totalObjects;
  unsigned char totalObstacles;
  unsigned char totalWaypoints;
  unsigned char verbosity;
  double visabilityModifier;
  double warningDistance;
  emxArray_char_T *outputFile;
} struct0_T;

#endif                                 /*typedef_struct0_T*/
#endif

/*
 * File trailer for OMAS_updateGlobalStates_types.h
 *
 * [EOF]
 */
