/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_getTrajectoryData_emxAPI.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Nov-2018 17:59:35
 */

/* Include Files */
#include <stdlib.h>
#include "rt_nonfinite.h"
#include "OMAS_getTrajectoryData.h"
#include "OMAS_getTrajectoryData_emxAPI.h"
#include "OMAS_getTrajectoryData_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : int numDimensions
 *                int *size
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T1(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : int numDimensions
 *                int *size
 * Return Type  : emxArray_uint8_T *
 */
emxArray_uint8_T *emxCreateND_uint8_T(int numDimensions, int *size)
{
  emxArray_uint8_T *emx;
  int numEl;
  int i;
  emxInit_uint8_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (unsigned char *)calloc((unsigned int)numEl, sizeof(unsigned char));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : double *data
 *                int numDimensions
 *                int *size
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateWrapperND_real_T(double *data, int numDimensions, int *
  size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T1(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : unsigned char *data
 *                int numDimensions
 *                int *size
 * Return Type  : emxArray_uint8_T *
 */
emxArray_uint8_T *emxCreateWrapperND_uint8_T(unsigned char *data, int
  numDimensions, int *size)
{
  emxArray_uint8_T *emx;
  int numEl;
  int i;
  emxInit_uint8_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : double *data
 *                int rows
 *                int cols
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T1(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : unsigned char *data
 *                int rows
 *                int cols
 * Return Type  : emxArray_uint8_T *
 */
emxArray_uint8_T *emxCreateWrapper_uint8_T(unsigned char *data, int rows, int
  cols)
{
  emxArray_uint8_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_uint8_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : int rows
 *                int cols
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreate_real_T(int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T1(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : int rows
 *                int cols
 * Return Type  : emxArray_uint8_T *
 */
emxArray_uint8_T *emxCreate_uint8_T(int rows, int cols)
{
  emxArray_uint8_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_uint8_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (unsigned char *)calloc((unsigned int)numEl, sizeof(unsigned char));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : emxArray_real_T *emxArray
 * Return Type  : void
 */
void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}

/*
 * Arguments    : emxArray_uint8_T *emxArray
 * Return Type  : void
 */
void emxDestroyArray_uint8_T(emxArray_uint8_T *emxArray)
{
  emxFree_uint8_T(&emxArray);
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T1(pEmxArray, numDimensions);
}

/*
 * Arguments    : emxArray_uint8_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInitArray_uint8_T(emxArray_uint8_T **pEmxArray, int numDimensions)
{
  emxInit_uint8_T(pEmxArray, numDimensions);
}

/*
 * File trailer for OMAS_getTrajectoryData_emxAPI.c
 *
 * [EOF]
 */
