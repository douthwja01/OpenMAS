/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_updateGlobalStates_emxutil.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

/* Include Files */
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "OMAS_updateGlobalStates.h"
#include "OMAS_updateGlobalStates_emxutil.h"

/* Function Declarations */
static void emxCopyMatrix_real32_T(float dst[3], const float src[3]);
static void emxCopyMatrix_real_T(double dst[10], const double src[10]);
static void emxCopyMatrix_real_T1(double dst[9], const double src[9]);
static void emxCopy_boolean_T(emxArray_boolean_T **dst, emxArray_boolean_T *
  const *src);
static void emxCopy_char_T(emxArray_char_T **dst, emxArray_char_T * const *src);
static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src);
static void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
static void emxFreeStruct_struct2_T(struct2_T *pStruct);
static void emxInitStruct_struct2_T(struct2_T *pStruct);

/* Function Definitions */

/*
 * Arguments    : float dst[3]
 *                const float src[3]
 * Return Type  : void
 */
static void emxCopyMatrix_real32_T(float dst[3], const float src[3])
{
  int i;
  for (i = 0; i < 3; i++) {
    dst[i] = src[i];
  }
}

/*
 * Arguments    : double dst[10]
 *                const double src[10]
 * Return Type  : void
 */
static void emxCopyMatrix_real_T(double dst[10], const double src[10])
{
  int i;
  for (i = 0; i < 10; i++) {
    dst[i] = src[i];
  }
}

/*
 * Arguments    : double dst[9]
 *                const double src[9]
 * Return Type  : void
 */
static void emxCopyMatrix_real_T1(double dst[9], const double src[9])
{
  int i;
  for (i = 0; i < 9; i++) {
    dst[i] = src[i];
  }
}

/*
 * Arguments    : emxArray_boolean_T **dst
 *                emxArray_boolean_T * const *src
 * Return Type  : void
 */
static void emxCopy_boolean_T(emxArray_boolean_T **dst, emxArray_boolean_T *
  const *src)
{
  int numElDst;
  int numElSrc;
  int i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity_boolean_T(*dst, numElDst);
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

/*
 * Arguments    : emxArray_char_T **dst
 *                emxArray_char_T * const *src
 * Return Type  : void
 */
static void emxCopy_char_T(emxArray_char_T **dst, emxArray_char_T * const *src)
{
  int numElDst;
  int numElSrc;
  int i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity_char_T(*dst, numElDst);
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

/*
 * Arguments    : emxArray_real_T **dst
 *                emxArray_real_T * const *src
 * Return Type  : void
 */
static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src)
{
  int numElDst;
  int numElSrc;
  int i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity_real_T(*dst, numElDst);
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

/*
 * Arguments    : emxArray_real_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
static void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, sizeof(double));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(double) * oldNumel);
      if (emxArray->canFreeData) {
        free((void *)emxArray->data);
      }
    }

    emxArray->data = (double *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : struct2_T *pStruct
 * Return Type  : void
 */
static void emxFreeStruct_struct2_T(struct2_T *pStruct)
{
  emxFree_real_T(&pStruct->timeVector);
}

/*
 * Arguments    : struct2_T *pStruct
 * Return Type  : void
 */
static void emxInitStruct_struct2_T(struct2_T *pStruct)
{
  emxInit_real_T(&pStruct->timeVector, 2);
}

/*
 * Arguments    : struct1_T *dst
 *                const struct1_T *src
 * Return Type  : void
 */
void emxCopyStruct_struct1_T(struct1_T *dst, const struct1_T *src)
{
  dst->objectID = src->objectID;
  emxCopy_char_T(&dst->name, &src->name);
  emxCopy_char_T(&dst->class, &src->class);
  dst->type = src->type;
  emxCopyMatrix_real32_T(dst->colour, src->colour);
  emxCopy_char_T(&dst->symbol, &src->symbol);
  dst->radius = src->radius;
  dst->detectionRadius = src->detectionRadius;
  dst->idleStatus = src->idleStatus;
  emxCopyMatrix_real_T(dst->globalState, src->globalState);
  emxCopyMatrix_real_T1(dst->R, src->R);
  emxCopy_real_T(&dst->relativePositions, &src->relativePositions);
  emxCopy_boolean_T(&dst->objectStatus, &src->objectStatus);
}

/*
 * Arguments    : emxArray_boolean_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, sizeof(boolean_T));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(boolean_T) * oldNumel);
      if (emxArray->canFreeData) {
        free((void *)emxArray->data);
      }
    }

    emxArray->data = (boolean_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : emxArray_char_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, sizeof(char));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(char) * oldNumel);
      if (emxArray->canFreeData) {
        free((void *)emxArray->data);
      }
    }

    emxArray->data = (char *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : struct0_T *pStruct
 * Return Type  : void
 */
void emxFreeStruct_struct0_T(struct0_T *pStruct)
{
  emxFree_struct1_T(&pStruct->OBJECTS);
  emxFreeStruct_struct2_T(&pStruct->TIME);
  emxFree_uint8_T(&pStruct->globalIDvector);
  emxFree_char_T(&pStruct->outputPath);
  emxFree_char_T(&pStruct->phase);
  emxFree_char_T(&pStruct->systemFile);
  emxFree_char_T(&pStruct->outputFile);
}

/*
 * Arguments    : struct1_T *pStruct
 * Return Type  : void
 */
void emxFreeStruct_struct1_T(struct1_T *pStruct)
{
  emxFree_char_T(&pStruct->name);
  emxFree_char_T(&pStruct->class);
  emxFree_char_T(&pStruct->symbol);
  emxFree_real_T(&pStruct->relativePositions);
  emxFree_boolean_T(&pStruct->objectStatus);
}

/*
 * Arguments    : emxArray_boolean_T **pEmxArray
 * Return Type  : void
 */
void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if (((*pEmxArray)->data != (boolean_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_char_T **pEmxArray
 * Return Type  : void
 */
void emxFree_char_T(emxArray_char_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_char_T *)NULL) {
    if (((*pEmxArray)->data != (char *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_char_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 * Return Type  : void
 */
void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_struct1_T **pEmxArray
 * Return Type  : void
 */
void emxFree_struct1_T(emxArray_struct1_T **pEmxArray)
{
  int numEl;
  int i;
  if (*pEmxArray != (emxArray_struct1_T *)NULL) {
    if ((*pEmxArray)->data != (struct1_T *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }

      for (i = 0; i < numEl; i++) {
        emxFreeStruct_struct1_T(&(*pEmxArray)->data[i]);
      }

      if ((*pEmxArray)->canFreeData) {
        free((void *)(*pEmxArray)->data);
      }
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_struct1_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_uint8_T **pEmxArray
 * Return Type  : void
 */
void emxFree_uint8_T(emxArray_uint8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_uint8_T *)NULL) {
    if (((*pEmxArray)->data != (unsigned char *)NULL) && (*pEmxArray)
        ->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_uint8_T *)NULL;
  }
}

/*
 * Arguments    : struct0_T *pStruct
 * Return Type  : void
 */
void emxInitStruct_struct0_T(struct0_T *pStruct)
{
  emxInit_struct1_T1(&pStruct->OBJECTS, 2);
  emxInitStruct_struct2_T(&pStruct->TIME);
  emxInit_uint8_T(&pStruct->globalIDvector, 2);
  emxInit_char_T(&pStruct->outputPath, 2);
  emxInit_char_T(&pStruct->phase, 2);
  emxInit_char_T(&pStruct->systemFile, 2);
  emxInit_char_T(&pStruct->outputFile, 2);
}

/*
 * Arguments    : struct1_T *pStruct
 * Return Type  : void
 */
void emxInitStruct_struct1_T(struct1_T *pStruct)
{
  emxInit_char_T(&pStruct->name, 2);
  emxInit_char_T(&pStruct->class, 2);
  emxInit_char_T(&pStruct->symbol, 2);
  emxInit_real_T(&pStruct->relativePositions, 2);
  emxInit_boolean_T(&pStruct->objectStatus, 2);
}

/*
 * Arguments    : emxArray_boolean_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions)
{
  emxArray_boolean_T *emxArray;
  int i;
  *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
  emxArray = *pEmxArray;
  emxArray->data = (boolean_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_char_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxArray_char_T *emxArray;
  int i;
  *pEmxArray = (emxArray_char_T *)malloc(sizeof(emxArray_char_T));
  emxArray = *pEmxArray;
  emxArray->data = (char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_struct1_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_struct1_T1(emxArray_struct1_T **pEmxArray, int numDimensions)
{
  emxArray_struct1_T *emxArray;
  int i;
  *pEmxArray = (emxArray_struct1_T *)malloc(sizeof(emxArray_struct1_T));
  emxArray = *pEmxArray;
  emxArray->data = (struct1_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_uint8_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int numDimensions)
{
  emxArray_uint8_T *emxArray;
  int i;
  *pEmxArray = (emxArray_uint8_T *)malloc(sizeof(emxArray_uint8_T));
  emxArray = *pEmxArray;
  emxArray->data = (unsigned char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * File trailer for OMAS_updateGlobalStates_emxutil.c
 *
 * [EOF]
 */
