/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_updateGlobalStates.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

#ifndef OMAS_UPDATEGLOBALSTATES_H
#define OMAS_UPDATEGLOBALSTATES_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "OMAS_updateGlobalStates_types.h"

/* Function Declarations */
extern void OMAS_updateGlobalStates(const struct0_T *SIM, unsigned char objectID,
  const double globalVelocity_k[3], const double quaternion_k[4], boolean_T
  idleStatus_k, struct1_T *METAObjUpdate);

#endif

/*
 * File trailer for OMAS_updateGlobalStates.h
 *
 * [EOF]
 */
