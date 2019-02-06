/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_updateGlobalStates.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "OMAS_updateGlobalStates.h"
#include "OMAS_geometry.h"
#include "OMAS_updateGlobalStates_emxutil.h"
#include <stdio.h>

/* Function Definitions */

/*
 * This function reallocates the global properties of each entity to the
 *  META.OBJECT series to allow faster reference and increased independance
 *  of the main cycle from the object cycles.
 * Arguments    : const struct0_T *SIM
 *                unsigned char objectID
 *                const double globalVelocity_k[3]
 *                const double quaternion_k[4]
 *                boolean_T idleStatus_k
 *                struct1_T *METAObjUpdate
 * Return Type  : void
 */
void OMAS_updateGlobalStates(const struct0_T *SIM, unsigned char objectID, const
  double globalVelocity_k[3], const double quaternion_k[4], boolean_T
  idleStatus_k, struct1_T *METAObjUpdate)
{
  int i;
  boolean_T y;
  boolean_T b[4];
  boolean_T exitg1;
  emxArray_char_T *charStr;
  emxArray_boolean_T *logicalIndices;
  unsigned char u0;
  int i0;
  signed char varargin_1;
  double b_index;
  double b_quaternion_k[4];

  /* %%%%%%%%%%%%%%%%%%%%%%%% UPDATE PROCEDURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  UPDATE THE GLOBAL STATE VECTORS FROM ENTITY.VIRTUAL */
  /*  CONFIRM OBJECT VARIABLES */
  for (i = 0; i < 4; i++) {
    b[i] = rtIsNaN(quaternion_k[i]);
  }

  y = false;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i < 4)) {
    if (b[i]) {
      y = true;
      exitg1 = true;
    } else {
      i++;
    }
  }

  if (y) {
    emxInit_char_T(&charStr, 2);
    u0 = objectID;
    if (objectID > 127) {
      u0 = 127U;
    }

    varargin_1 = (signed char)u0;
    i = (int)snprintf(NULL, 0,
                      "Quaternion for objectID %d is invalid; q = [%f %f %f %f]",
                      varargin_1, quaternion_k[0], quaternion_k[1],
                      quaternion_k[2], quaternion_k[3]) + 1;
    i0 = charStr->size[0] * charStr->size[1];
    charStr->size[0] = 1;
    charStr->size[1] = i;
    emxEnsureCapacity_char_T(charStr, i0);
    snprintf(&charStr->data[0], (size_t)i,
             "Quaternion for objectID %d is invalid; q = [%f %f %f %f]",
             varargin_1, quaternion_k[0], quaternion_k[1], quaternion_k[2],
             quaternion_k[3]);
    emxFree_char_T(&charStr);
  }

  emxInit_boolean_T(&logicalIndices, 2);

  /*  IDENTIFY THE ASSOCIATED META OBJECT  */
  i0 = logicalIndices->size[0] * logicalIndices->size[1];
  logicalIndices->size[0] = 1;
  logicalIndices->size[1] = SIM->globalIDvector->size[1];
  emxEnsureCapacity_boolean_T(logicalIndices, i0);
  i = SIM->globalIDvector->size[0] * SIM->globalIDvector->size[1];
  for (i0 = 0; i0 < i; i0++) {
    logicalIndices->data[i0] = (SIM->globalIDvector->data[i0] == objectID);
  }

  /*  The logical position of the ID */
  b_index = rtInf;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i <= logicalIndices->size[1] - 1)) {
    if (logicalIndices->data[i]) {
      b_index = 1.0 + (double)i;
      exitg1 = true;
    } else {
      i++;
    }
  }

  emxFree_boolean_T(&logicalIndices);

  /*  EXTRACT THE META STRUCTURE TO BE UPDATED  */
  emxCopyStruct_struct1_T(METAObjUpdate, &SIM->OBJECTS->data[SIM->OBJECTS->size
    [0] * ((int)b_index - 1)]);

  /*  Get the current META object associated with 'entity' */
  /*  //////////////////////// UPDATE META.OBJECT ///////////////////////////// */
  /*  UPDATE THE META ROTATION DEFINING FIXED LOCAL >> ROTATED IN THE GLOBAL */
  for (i = 0; i < 4; i++) {
    b_quaternion_k[i] = quaternion_k[i];
  }

  c_OMAS_geometry_quaternionToRot(b_quaternion_k, METAObjUpdate->R);

  /*  METAObjUpdate.R = quat2rotm(quaternion_k');   */
  /*  UPDATE THE GLOBAL POSITION */
  /*  Calculate the new global position */
  /*  REBUILD GLOBAL STATES (ENTITY & META) */
  for (i0 = 0; i0 < 3; i0++) {
    METAObjUpdate->globalState[i0] = SIM->OBJECTS->data[SIM->OBJECTS->size[0] *
      ((int)b_index - 1)].globalState[i0] + globalVelocity_k[i0] * SIM->TIME.dt;
    METAObjUpdate->globalState[i0 + 3] = globalVelocity_k[i0];
  }

  for (i0 = 0; i0 < 4; i0++) {
    METAObjUpdate->globalState[i0 + 6] = quaternion_k[i0];
  }

  /*  DETERMINE IF ENTITY HAS INDICATED THAT TASK IS COMPLETE */
  METAObjUpdate->idleStatus = idleStatus_k;
}

/*
 * File trailer for OMAS_updateGlobalStates.c
 *
 * [EOF]
 */
