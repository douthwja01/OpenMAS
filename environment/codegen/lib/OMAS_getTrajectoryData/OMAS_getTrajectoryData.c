/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_getTrajectoryData.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Nov-2018 17:59:35
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "OMAS_getTrajectoryData.h"
#include "OMAS_getTrajectoryData_emxutil.h"

/* Function Declarations */
static int div_nde_s32_floor(int numerator, int denominator);

/* Function Definitions */

/*
 * Arguments    : int numerator
 *                int denominator
 * Return Type  : int
 */
static int div_nde_s32_floor(int numerator, int denominator)
{
  int b_numerator;
  if (((numerator < 0) != (denominator < 0)) && (numerator % denominator != 0))
  {
    b_numerator = -1;
  } else {
    b_numerator = 0;
  }

  return numerator / denominator + b_numerator;
}

/*
 * This is the generic function for returning the state trajectory data for
 *  a given agent ID number.
 *  INPUTS:
 *  trajectoryMatrix - a complete matrix of object global states [objects*n x steps].
 *  objectNumber     - The the index of desired objects ID in the SIM.OBJECTS list.
 *  step             - The requested simulation step.
 *   + >steps        - Return states for all timesteps.
 *   + nan           - Return states for all timesteps until the agent way idle.
 *   + step          - Return a specific state at a given timestep.
 * Arguments    : const emxArray_real_T *trajectoryMatrix
 *                const emxArray_uint8_T *globalIDvector
 *                unsigned char objectID
 *                float step
 *                emxArray_real_T *objectStateData
 *                unsigned int *indexOfLastState
 * Return Type  : void
 */
void OMAS_getTrajectoryData(const emxArray_real_T *trajectoryMatrix, const
  emxArray_uint8_T *globalIDvector, unsigned char objectID, float step,
  emxArray_real_T *objectStateData, unsigned int *indexOfLastState)
{
  emxArray_boolean_T *logicalIndices;
  int i0;
  int i1;
  double indexPosition;
  int i;
  boolean_T exitg1;
  emxArray_uint16_T *indexSet;
  double y;
  emxArray_int32_T *b_y;
  int b_i1;
  unsigned short u0;
  emxArray_boolean_T *logicalMatrix;
  emxArray_real_T *b_trajectoryMatrix;
  emxArray_real_T *b_objectStateData;
  int i2;
  int iy;
  int b_i2;
  unsigned int outsize[2];
  boolean_T c_y;
  boolean_T b0;
  int b_indexOfLastState;
  emxInit_boolean_T(&logicalIndices, 2);

  /*  OPENMAS TRAJECTORY DATA EXTRACTION (simulation_getTrajectory.m) %%%%%%%% */
  /*  This function is designed to provide a utility for extracting object */
  /*  trajectory data from the raw DATA.globalTrajectories data set. */
  /*  Author: James A. Douthwaite */
  /*  Fetch states from Trajectory matrix */
  /*  DETERMINE THE NUMBER OF OUTPUT STATES */
  /*  Determing the number of system states */
  /*  INPUT HANDLING */
  /*  DEFINE THE INDEX OF THE OBJECT FROM ITS OBJECT ID  */
  i0 = logicalIndices->size[0] * logicalIndices->size[1];
  logicalIndices->size[0] = 1;
  logicalIndices->size[1] = globalIDvector->size[1];
  emxEnsureCapacity_boolean_T(logicalIndices, i0);
  i1 = globalIDvector->size[0] * globalIDvector->size[1];
  for (i0 = 0; i0 < i1; i0++) {
    logicalIndices->data[i0] = (globalIDvector->data[i0] == objectID);
  }

  /*  The logical position of the ID */
  indexPosition = rtInf;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i <= logicalIndices->size[1] - 1)) {
    if (logicalIndices->data[i]) {
      indexPosition = 1.0 + (double)i;
      exitg1 = true;
    } else {
      i++;
    }
  }

  emxInit_uint16_T(&indexSet, 2);

  /* INFER INDEX PROPERTIES */
  y = (double)trajectoryMatrix->size[0] / 10.0;
  i0 = indexSet->size[0] * indexSet->size[1];
  indexSet->size[0] = 2;
  indexSet->size[1] = (int)y;
  emxEnsureCapacity_uint16_T(indexSet, i0);
  i1 = (int)y << 1;
  for (i0 = 0; i0 < i1; i0++) {
    indexSet->data[i0] = 0U;
  }

  emxInit_int32_T(&b_y, 2);
  if (trajectoryMatrix->size[0] - 10 < 0) {
    i0 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = 0;
    emxEnsureCapacity_int32_T(b_y, i0);
  } else {
    i0 = trajectoryMatrix->size[0] - 10;
    b_i1 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = div_nde_s32_floor(i0, 10) + 1;
    emxEnsureCapacity_int32_T(b_y, b_i1);
    i1 = div_nde_s32_floor(i0, 10);
    for (i0 = 0; i0 <= i1; i0++) {
      b_y->data[b_y->size[0] * i0] = 10 * i0;
    }
  }

  i1 = b_y->size[1];
  for (i0 = 0; i0 < i1; i0++) {
    y = (double)b_y->data[b_y->size[0] * i0] + 1.0;
    if (y < 65536.0) {
      if (y >= 0.0) {
        u0 = (unsigned short)y;
      } else {
        u0 = 0U;
      }
    } else {
      u0 = MAX_uint16_T;
    }

    indexSet->data[indexSet->size[0] * i0] = u0;
  }

  /*  Declare the state indices */
  if (trajectoryMatrix->size[0] < 10) {
    i0 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = 0;
    emxEnsureCapacity_int32_T(b_y, i0);
  } else {
    i0 = trajectoryMatrix->size[0] - 10;
    b_i1 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = div_nde_s32_floor(i0, 10) + 1;
    emxEnsureCapacity_int32_T(b_y, b_i1);
    i1 = div_nde_s32_floor(i0, 10);
    for (i0 = 0; i0 <= i1; i0++) {
      b_y->data[b_y->size[0] * i0] = 10 + 10 * i0;
    }
  }

  i1 = b_y->size[1];
  for (i0 = 0; i0 < i1; i0++) {
    b_i1 = b_y->data[b_y->size[0] * i0];
    if (b_i1 < 0) {
      b_i1 = 0;
    } else {
      if (b_i1 > 65535) {
        b_i1 = 65535;
      }
    }

    indexSet->data[1 + indexSet->size[0] * i0] = (unsigned short)b_i1;
  }

  emxFree_int32_T(&b_y);

  /*  DEFINE THE STATE LIMITS (INDICES) THAT  */
  /*  RETURN ALL TIMESERIES DATA */
  if ((double)step > trajectoryMatrix->size[1]) {
    /*  FULL TIME SERIES */
    if (indexSet->data[indexSet->size[0] * ((int)indexPosition - 1)] >
        indexSet->data[1 + indexSet->size[0] * ((int)indexPosition - 1)]) {
      i0 = 0;
      b_i1 = 0;
    } else {
      i0 = indexSet->data[indexSet->size[0] * ((int)indexPosition - 1)] - 1;
      b_i1 = indexSet->data[1 + indexSet->size[0] * ((int)indexPosition - 1)];
    }

    i1 = trajectoryMatrix->size[1];
    i2 = objectStateData->size[0] * objectStateData->size[1];
    objectStateData->size[0] = b_i1 - i0;
    objectStateData->size[1] = i1;
    emxEnsureCapacity_real_T(objectStateData, i2);
    for (i2 = 0; i2 < i1; i2++) {
      iy = b_i1 - i0;
      for (b_i2 = 0; b_i2 < iy; b_i2++) {
        objectStateData->data[b_i2 + objectStateData->size[0] * i2] =
          trajectoryMatrix->data[(i0 + b_i2) + trajectoryMatrix->size[0] * i2];
      }
    }

    /*  Otherwise extract the complete timeseries */
    i0 = trajectoryMatrix->size[1];
    *indexOfLastState = (unsigned int)i0;

    /*  Indices of the last state */
  } else {
    /*  assert(isnumeric(step),'Step provided must be an integer number.'); */
    /*  RETURN ALL STATES UPTO THE POINT OF INACTIVITY (INDICATED BY NaNs) */
    emxInit_boolean_T(&logicalMatrix, 2);
    emxInit_real_T(&b_trajectoryMatrix, 1);
    emxInit_real_T1(&b_objectStateData, 2);
    if (rtIsNaNF(step)) {
      /*  RETURN ONLY VALID STATES (ALL STATES UPTO TERMINATION 'NaN') */
      if (indexSet->data[indexSet->size[0] * ((int)indexPosition - 1)] >
          indexSet->data[1 + indexSet->size[0] * ((int)indexPosition - 1)]) {
        i0 = 1;
        b_i1 = 1;
      } else {
        i0 = indexSet->data[indexSet->size[0] * ((int)indexPosition - 1)];
        b_i1 = indexSet->data[1 + indexSet->size[0] * ((int)indexPosition - 1)]
          + 1;
      }

      i1 = trajectoryMatrix->size[1];
      i2 = objectStateData->size[0] * objectStateData->size[1];
      objectStateData->size[0] = b_i1 - i0;
      objectStateData->size[1] = i1;
      emxEnsureCapacity_real_T(objectStateData, i2);
      for (i2 = 0; i2 < i1; i2++) {
        iy = b_i1 - i0;
        for (b_i2 = 0; b_i2 < iy; b_i2++) {
          objectStateData->data[b_i2 + objectStateData->size[0] * i2] =
            trajectoryMatrix->data[((i0 + b_i2) + trajectoryMatrix->size[0] * i2)
            - 1];
        }
      }

      /*  Full timeseries */
      /*  IF THE AGENT IS ACTIVE AT ALL TIMESTEPS */
      i1 = trajectoryMatrix->size[1];
      i2 = logicalMatrix->size[0] * logicalMatrix->size[1];
      logicalMatrix->size[0] = 3;
      logicalMatrix->size[1] = i1;
      emxEnsureCapacity_boolean_T(logicalMatrix, i2);
      for (i2 = 0; i2 < i1; i2++) {
        for (b_i2 = 0; b_i2 < 3; b_i2++) {
          logicalMatrix->data[b_i2 + logicalMatrix->size[0] * i2] = rtIsNaN
            (objectStateData->data[b_i2 + objectStateData->size[0] * i2]);
        }
      }

      for (i2 = 0; i2 < 2; i2++) {
        outsize[i2] = (unsigned int)logicalMatrix->size[i2];
      }

      i2 = logicalIndices->size[0] * logicalIndices->size[1];
      logicalIndices->size[0] = 1;
      logicalIndices->size[1] = (int)outsize[1];
      emxEnsureCapacity_boolean_T(logicalIndices, i2);
      i1 = (int)outsize[1];
      for (i2 = 0; i2 < i1; i2++) {
        logicalIndices->data[i2] = false;
      }

      b_i2 = 1;
      iy = -1;
      for (i = 1; i <= logicalMatrix->size[1]; i++) {
        i1 = b_i2;
        b_i2 += 3;
        iy++;
        exitg1 = false;
        while ((!exitg1) && (i1 <= b_i2 - 1)) {
          b0 = !logicalMatrix->data[i1 - 1];
          if (!b0) {
            logicalIndices->data[iy] = true;
            exitg1 = true;
          } else {
            i1++;
          }
        }
      }

      c_y = false;
      b_i2 = 1;
      exitg1 = false;
      while ((!exitg1) && (b_i2 <= logicalIndices->size[1])) {
        b0 = !logicalIndices->data[b_i2 - 1];
        if (!b0) {
          c_y = true;
          exitg1 = true;
        } else {
          b_i2++;
        }
      }

      if (!c_y) {
        i0 = trajectoryMatrix->size[1];
        *indexOfLastState = (unsigned int)i0;
      } else {
        /*  THERE ARE NAN OCCURANCES */
        b_indexOfLastState = logicalMatrix->size[1];
        i2 = logicalMatrix->size[1];
        i = 0;
        exitg1 = false;
        while ((!exitg1) && (i <= 1 - i2)) {
          b_i2 = i2 + i;
          y = logicalMatrix->data[logicalMatrix->size[0] * (b_i2 - 1)];
          for (iy = 0; iy < 2; iy++) {
            y += (double)logicalMatrix->data[(iy + logicalMatrix->size[0] *
              (b_i2 - 1)) + 1];
          }

          if (y != 3.0) {
            b_indexOfLastState = logicalMatrix->size[1] - b_i2;
            exitg1 = true;
          } else {
            i++;
          }
        }

        if (1 > b_indexOfLastState) {
          i1 = 0;
        } else {
          i1 = b_indexOfLastState;
        }

        b_i2 = b_i1 - i0;
        i0 = b_objectStateData->size[0] * b_objectStateData->size[1];
        b_objectStateData->size[0] = b_i2;
        b_objectStateData->size[1] = i1;
        emxEnsureCapacity_real_T(b_objectStateData, i0);
        for (i0 = 0; i0 < i1; i0++) {
          for (b_i1 = 0; b_i1 < b_i2; b_i1++) {
            b_objectStateData->data[b_i1 + b_objectStateData->size[0] * i0] =
              objectStateData->data[b_i1 + objectStateData->size[0] * i0];
          }
        }

        i0 = objectStateData->size[0] * objectStateData->size[1];
        objectStateData->size[0] = b_objectStateData->size[0];
        objectStateData->size[1] = b_objectStateData->size[1];
        emxEnsureCapacity_real_T(objectStateData, i0);
        i1 = b_objectStateData->size[1];
        for (i0 = 0; i0 < i1; i0++) {
          iy = b_objectStateData->size[0];
          for (b_i1 = 0; b_i1 < iy; b_i1++) {
            objectStateData->data[b_i1 + objectStateData->size[0] * i0] =
              b_objectStateData->data[b_i1 + b_objectStateData->size[0] * i0];
          }
        }

        /*  Final valid object state */
        if (b_indexOfLastState < 0) {
          b_indexOfLastState = 0;
        }

        *indexOfLastState = (unsigned int)b_indexOfLastState;
      }
    } else {
      /*  SPECIFIC STEP */
      if (indexSet->data[indexSet->size[0] * ((int)indexPosition - 1)] >
          indexSet->data[1 + indexSet->size[0] * ((int)indexPosition - 1)]) {
        i0 = 1;
        b_i1 = 1;
      } else {
        i0 = indexSet->data[indexSet->size[0] * ((int)indexPosition - 1)];
        b_i1 = indexSet->data[1 + indexSet->size[0] * ((int)indexPosition - 1)]
          + 1;
      }

      i2 = b_trajectoryMatrix->size[0];
      b_trajectoryMatrix->size[0] = b_i1 - i0;
      emxEnsureCapacity_real_T1(b_trajectoryMatrix, i2);
      i1 = b_i1 - i0;
      for (i2 = 0; i2 < i1; i2++) {
        b_trajectoryMatrix->data[i2] = trajectoryMatrix->data[((i0 + i2) +
          trajectoryMatrix->size[0] * ((int)step - 1)) - 1];
      }

      b_i2 = b_i1 - i0;
      i0 = objectStateData->size[0] * objectStateData->size[1];
      objectStateData->size[0] = b_i2;
      objectStateData->size[1] = 1;
      emxEnsureCapacity_real_T(objectStateData, i0);
      for (i0 = 0; i0 < 1; i0++) {
        for (b_i1 = 0; b_i1 < b_i2; b_i1++) {
          objectStateData->data[b_i1] = b_trajectoryMatrix->data[b_i1];
        }
      }

      /*  A specific time step has been specified */
      *indexOfLastState = (unsigned int)step;
    }

    emxFree_real_T(&b_objectStateData);
    emxFree_real_T(&b_trajectoryMatrix);
    emxFree_boolean_T(&logicalMatrix);
  }

  emxFree_uint16_T(&indexSet);
  emxFree_boolean_T(&logicalIndices);

  /*  if ischar(step)  */
  /*      if strncmp(step,'last',4) */
  /*          % LAST STATE REQUESTED */
  /*          objectStateData = trajectoryMatrix(objectStateIndices(1):objectStateIndices(2),:);   % Full timeseries */
  /*          validStates = ~any(isnan(objectStateData(1:3,:)));                            % non-NaN states */
  /*          indicesOfFinalState = find(validStates,1,'last'); */
  /*          objectStateData = objectStateData(:,indicesOfFinalState);                              % Final valid object state */
  /*      elseif strncmp(step,'valid',5) */
}

/*
 * File trailer for OMAS_getTrajectoryData.c
 *
 * [EOF]
 */
