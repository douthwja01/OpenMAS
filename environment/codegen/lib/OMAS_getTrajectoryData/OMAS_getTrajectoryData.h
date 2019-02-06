/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_getTrajectoryData.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Nov-2018 17:59:35
 */

#ifndef OMAS_GETTRAJECTORYDATA_H
#define OMAS_GETTRAJECTORYDATA_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "OMAS_getTrajectoryData_types.h"

/* Function Declarations */
extern void OMAS_getTrajectoryData(const emxArray_real_T *trajectoryMatrix,
  const emxArray_uint8_T *globalIDvector, unsigned char objectID, float step,
  emxArray_real_T *objectStateData, unsigned int *indexOfLastState);

#endif

/*
 * File trailer for OMAS_getTrajectoryData.h
 *
 * [EOF]
 */
