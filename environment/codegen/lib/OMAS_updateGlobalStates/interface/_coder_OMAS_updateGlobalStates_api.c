/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_OMAS_updateGlobalStates_api.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

/* Include Files */
#include <string.h>
#include "tmwtypes.h"
#include "_coder_OMAS_updateGlobalStates_api.h"
#include "_coder_OMAS_updateGlobalStates_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131466U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "OMAS_updateGlobalStates",           /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static boolean_T ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static void bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[10]);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_struct1_T *y);
static void cb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[9]);
static uint8_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void db_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_char_T *y);
static void eb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_boolean_T *ret);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *SIM, const
  char_T *identifier, struct0_T *y);
static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const struct1_T *u);
static void emxEnsureCapacity_boolean_T(const emlrtStack *sp, emxArray_boolean_T
  *emxArray, int32_T oldNumel);
static void emxEnsureCapacity_char_T(const emlrtStack *sp, emxArray_char_T
  *emxArray, int32_T oldNumel);
static void emxEnsureCapacity_real_T(const emlrtStack *sp, emxArray_real_T
  *emxArray, int32_T oldNumel);
static void emxEnsureCapacity_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  *emxArray, int32_T oldNumel);
static void emxEnsureCapacity_uint8_T(const emlrtStack *sp, emxArray_uint8_T
  *emxArray, int32_T oldNumel);
static void emxExpand_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  *emxArray, int32_T fromIndex, int32_T toIndex);
static void emxFreeStruct_struct0_T(const emlrtStack *sp, struct0_T *pStruct);
static void emxFreeStruct_struct1_T(const emlrtStack *sp, struct1_T *pStruct);
static void emxFreeStruct_struct2_T(const emlrtStack *sp, struct2_T *pStruct);
static void emxFree_boolean_T(const emlrtStack *sp, emxArray_boolean_T
  **pEmxArray);
static void emxFree_char_T(const emlrtStack *sp, emxArray_char_T **pEmxArray);
static void emxFree_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray);
static void emxFree_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  **pEmxArray);
static void emxFree_uint8_T(const emlrtStack *sp, emxArray_uint8_T **pEmxArray);
static void emxInitStruct_struct0_T(const emlrtStack *sp, struct0_T *pStruct,
  boolean_T doPush);
static void emxInitStruct_struct1_T(const emlrtStack *sp, struct1_T *pStruct,
  boolean_T doPush);
static void emxInitStruct_struct2_T(const emlrtStack *sp, struct2_T *pStruct,
  boolean_T doPush);
static void emxInit_boolean_T(const emlrtStack *sp, emxArray_boolean_T
  **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void emxInit_char_T(const emlrtStack *sp, emxArray_char_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static void emxInit_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void emxInit_uint8_T(const emlrtStack *sp, emxArray_uint8_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static void emxTrim_struct1_T(const emlrtStack *sp, emxArray_struct1_T *emxArray,
  int32_T fromIndex, int32_T toIndex);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[3]);
static void fb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void gb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_uint8_T *ret);
static boolean_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId);
static real_T (*hb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[3];
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[10]);
static real_T (*ib_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4];
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[9]);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_boolean_T *y);
static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct2_T *y);
static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_uint8_T *y);
static uint8_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *objectID,
  const char_T *identifier);
static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *globalVelocity_k, const char_T *identifier))[3];
static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[3];
static real_T (*s_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *quaternion_k, const char_T *identifier))[4];
static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4];
static boolean_T u_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *idleStatus_k, const char_T *identifier);
static uint8_T v_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static void w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_char_T *ret);
static void x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[3]);
static real_T y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : boolean_T
 */
static boolean_T ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  boolean_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "logical", false, 0U, &dims);
  ret = *emlrtMxGetLogicals(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                struct0_T *y
 * Return Type  : void
 */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[19] = { "OBJECTS", "TIME", "conditionTolerance",
    "globalIDvector", "gui", "monteCarloMode", "outputPath", "phase",
    "systemFile", "threadPool", "totalAgents", "totalMiscs", "totalObjects",
    "totalObstacles", "totalWaypoints", "verbosity", "visabilityModifier",
    "warningDistance", "outputFile" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 19, fieldNames, 0U, &dims);
  thisId.fIdentifier = "OBJECTS";
  c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 0, "OBJECTS")),
                     &thisId, y->OBJECTS);
  thisId.fIdentifier = "TIME";
  m_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 1, "TIME")),
                     &thisId, &y->TIME);
  thisId.fIdentifier = "conditionTolerance";
  y->conditionTolerance = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    (sp, u, 0, 2, "conditionTolerance")), &thisId);
  thisId.fIdentifier = "globalIDvector";
  o_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 3,
    "globalIDvector")), &thisId, y->globalIDvector);
  thisId.fIdentifier = "gui";
  y->gui = h_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4,
    "gui")), &thisId);
  thisId.fIdentifier = "monteCarloMode";
  y->monteCarloMode = h_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 5, "monteCarloMode")), &thisId);
  thisId.fIdentifier = "outputPath";
  e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 6,
    "outputPath")), &thisId, y->outputPath);
  thisId.fIdentifier = "phase";
  e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 7, "phase")),
                     &thisId, y->phase);
  thisId.fIdentifier = "systemFile";
  e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 8,
    "systemFile")), &thisId, y->systemFile);
  thisId.fIdentifier = "threadPool";
  y->threadPool = h_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    9, "threadPool")), &thisId);
  thisId.fIdentifier = "totalAgents";
  y->totalAgents = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 10, "totalAgents")), &thisId);
  thisId.fIdentifier = "totalMiscs";
  y->totalMiscs = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    11, "totalMiscs")), &thisId);
  thisId.fIdentifier = "totalObjects";
  y->totalObjects = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 12, "totalObjects")), &thisId);
  thisId.fIdentifier = "totalObstacles";
  y->totalObstacles = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 13, "totalObstacles")), &thisId);
  thisId.fIdentifier = "totalWaypoints";
  y->totalWaypoints = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 14, "totalWaypoints")), &thisId);
  thisId.fIdentifier = "verbosity";
  y->verbosity = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    15, "verbosity")), &thisId);
  thisId.fIdentifier = "visabilityModifier";
  y->visabilityModifier = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    (sp, u, 0, 16, "visabilityModifier")), &thisId);
  thisId.fIdentifier = "warningDistance";
  y->warningDistance = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 17, "warningDistance")), &thisId);
  thisId.fIdentifier = "outputFile";
  e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 18,
    "outputFile")), &thisId, y->outputFile);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[10]
 * Return Type  : void
 */
static void bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[10])
{
  static const int32_T dims[1] = { 10 };

  int32_T i4;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  for (i4 = 0; i4 < 10; i4++) {
    ret[i4] = (*(real_T (*)[10])emlrtMxGetData(src))[i4];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                emxArray_struct1_T *y
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_struct1_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[13] = { "objectID", "name", "class", "type",
    "colour", "symbol", "radius", "detectionRadius", "idleStatus", "globalState",
    "R", "relativePositions", "objectStatus" };

  static const int32_T dims[2] = { 1, -1 };

  const boolean_T bv0[2] = { false, true };

  int32_T sizes[2];
  int32_T i0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckVsStructR2012b(sp, parentId, u, 13, fieldNames, 2U, dims, &bv0[0],
    sizes);
  i0 = y->size[0] * y->size[1];
  y->size[0] = sizes[0];
  y->size[1] = sizes[1];
  emxEnsureCapacity_struct1_T(sp, y, i0);
  for (i0 = 0; i0 < sizes[1]; i0++) {
    thisId.fIdentifier = "objectID";
    y->data[i0].objectID = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
      (sp, u, i0, 0, "objectID")), &thisId);
    thisId.fIdentifier = "name";
    e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 1, "name")),
                       &thisId, y->data[i0].name);
    thisId.fIdentifier = "class";
    e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 2, "class")),
                       &thisId, y->data[i0].class);
    thisId.fIdentifier = "type";
    y->data[i0].type = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
      u, i0, 3, "type")), &thisId);
    thisId.fIdentifier = "colour";
    f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 4, "colour")),
                       &thisId, y->data[i0].colour);
    thisId.fIdentifier = "symbol";
    e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 5, "symbol")),
                       &thisId, y->data[i0].symbol);
    thisId.fIdentifier = "radius";
    y->data[i0].radius = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
      (sp, u, i0, 6, "radius")), &thisId);
    thisId.fIdentifier = "detectionRadius";
    y->data[i0].detectionRadius = g_emlrt_marshallIn(sp, emlrtAlias
      (emlrtGetFieldR2017b(sp, u, i0, 7, "detectionRadius")), &thisId);
    thisId.fIdentifier = "idleStatus";
    y->data[i0].idleStatus = h_emlrt_marshallIn(sp, emlrtAlias
      (emlrtGetFieldR2017b(sp, u, i0, 8, "idleStatus")), &thisId);
    thisId.fIdentifier = "globalState";
    i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 9,
      "globalState")), &thisId, y->data[i0].globalState);
    thisId.fIdentifier = "R";
    j_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 10, "R")),
                       &thisId, y->data[i0].R);
    thisId.fIdentifier = "relativePositions";
    k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 11,
      "relativePositions")), &thisId, y->data[i0].relativePositions);
    thisId.fIdentifier = "objectStatus";
    l_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i0, 12,
      "objectStatus")), &thisId, y->data[i0].objectStatus);
  }

  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[9]
 * Return Type  : void
 */
static void cb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[9])
{
  static const int32_T dims[2] = { 3, 3 };

  int32_T i5;
  int32_T i6;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  for (i5 = 0; i5 < 3; i5++) {
    for (i6 = 0; i6 < 3; i6++) {
      ret[i6 + 3 * i5] = (*(real_T (*)[9])emlrtMxGetData(src))[i6 + 3 * i5];
    }
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : uint8_T
 */
static uint8_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  uint8_T y;
  y = v_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                emxArray_real_T *ret
 * Return Type  : void
 */
static void db_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { -1, 3 };

  const boolean_T bv2[2] = { true, false };

  int32_T iv4[2];
  int32_T i7;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv2[0],
    iv4);
  i7 = ret->size[0] * ret->size[1];
  ret->size[0] = iv4[0];
  ret->size[1] = iv4[1];
  emxEnsureCapacity_real_T(sp, ret, i7);
  emlrtImportArrayR2015b(sp, src, ret->data, 8, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                emxArray_char_T *y
 * Return Type  : void
 */
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_char_T *y)
{
  w_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                emxArray_boolean_T *ret
 * Return Type  : void
 */
static void eb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_boolean_T *ret)
{
  static const int32_T dims[2] = { -1, 4 };

  const boolean_T bv3[2] = { true, false };

  int32_T iv5[2];
  int32_T i8;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "logical", false, 2U, dims, &bv3[0],
    iv5);
  i8 = ret->size[0] * ret->size[1];
  ret->size[0] = iv5[0];
  ret->size[1] = iv5[1];
  emxEnsureCapacity_boolean_T(sp, ret, i8);
  emlrtImportArrayR2015b(sp, src, ret->data, 1, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *SIM
 *                const char_T *identifier
 *                struct0_T *y
 * Return Type  : void
 */
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *SIM, const
  char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(SIM), &thisId, y);
  emlrtDestroyArray(&SIM);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const struct1_T *u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const struct1_T *u)
{
  const mxArray *y;
  emxArray_char_T *b_u;
  static const char * sv0[13] = { "objectID", "name", "class", "type", "colour",
    "symbol", "radius", "detectionRadius", "idleStatus", "globalState", "R",
    "relativePositions", "objectStatus" };

  const mxArray *b_y;
  const mxArray *m0;
  int32_T i1;
  int32_T loop_ub;
  static const int32_T iv0[2] = { 1, 3 };

  real32_T *pData;
  static const int32_T iv1[1] = { 10 };

  real_T *b_pData;
  static const int32_T iv2[2] = { 3, 3 };

  emxArray_real_T *c_u;
  int32_T i;
  emxArray_boolean_T *d_u;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_char_T(sp, &b_u, 2, true);
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 13, sv0));
  b_y = NULL;
  m0 = emlrtCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
  *(uint8_T *)emlrtMxGetData(m0) = u->objectID;
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "objectID", b_y, 0);
  i1 = b_u->size[0] * b_u->size[1];
  b_u->size[0] = 1;
  b_u->size[1] = u->name->size[1];
  emxEnsureCapacity_char_T(sp, b_u, i1);
  loop_ub = u->name->size[0] * u->name->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    b_u->data[i1] = u->name->data[i1];
  }

  b_y = NULL;
  m0 = emlrtCreateCharArray(2, *(int32_T (*)[2])b_u->size);
  emlrtInitCharArrayR2013a(sp, u->name->size[1], m0, &b_u->data[0]);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "name", b_y, 1);
  i1 = b_u->size[0] * b_u->size[1];
  b_u->size[0] = 1;
  b_u->size[1] = u->class->size[1];
  emxEnsureCapacity_char_T(sp, b_u, i1);
  loop_ub = u->class->size[0] * u->class->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    b_u->data[i1] = u->class->data[i1];
  }

  b_y = NULL;
  m0 = emlrtCreateCharArray(2, *(int32_T (*)[2])b_u->size);
  emlrtInitCharArrayR2013a(sp, u->class->size[1], m0, &b_u->data[0]);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "class", b_y, 2);
  b_y = NULL;
  m0 = emlrtCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
  *(uint8_T *)emlrtMxGetData(m0) = u->type;
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "type", b_y, 3);
  b_y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m0);
  i1 = 0;
  for (loop_ub = 0; loop_ub < 3; loop_ub++) {
    pData[i1] = u->colour[loop_ub];
    i1++;
  }

  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "colour", b_y, 4);
  i1 = b_u->size[0] * b_u->size[1];
  b_u->size[0] = 1;
  b_u->size[1] = u->symbol->size[1];
  emxEnsureCapacity_char_T(sp, b_u, i1);
  loop_ub = u->symbol->size[0] * u->symbol->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    b_u->data[i1] = u->symbol->data[i1];
  }

  b_y = NULL;
  m0 = emlrtCreateCharArray(2, *(int32_T (*)[2])b_u->size);
  emlrtInitCharArrayR2013a(sp, u->symbol->size[1], m0, &b_u->data[0]);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "symbol", b_y, 5);
  b_y = NULL;
  m0 = emlrtCreateDoubleScalar(u->radius);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "radius", b_y, 6);
  b_y = NULL;
  m0 = emlrtCreateDoubleScalar(u->detectionRadius);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "detectionRadius", b_y, 7);
  b_y = NULL;
  m0 = emlrtCreateLogicalScalar(u->idleStatus);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "idleStatus", b_y, 8);
  b_y = NULL;
  m0 = emlrtCreateNumericArray(1, iv1, mxDOUBLE_CLASS, mxREAL);
  b_pData = emlrtMxGetPr(m0);
  i1 = 0;
  emxFree_char_T(sp, &b_u);
  for (loop_ub = 0; loop_ub < 10; loop_ub++) {
    b_pData[i1] = u->globalState[loop_ub];
    i1++;
  }

  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "globalState", b_y, 9);
  b_y = NULL;
  m0 = emlrtCreateNumericArray(2, iv2, mxDOUBLE_CLASS, mxREAL);
  b_pData = emlrtMxGetPr(m0);
  i1 = 0;
  for (loop_ub = 0; loop_ub < 3; loop_ub++) {
    for (i = 0; i < 3; i++) {
      b_pData[i1] = u->R[i + 3 * loop_ub];
      i1++;
    }
  }

  emxInit_real_T(sp, &c_u, 2, true);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "R", b_y, 10);
  i1 = c_u->size[0] * c_u->size[1];
  c_u->size[0] = u->relativePositions->size[0];
  c_u->size[1] = 3;
  emxEnsureCapacity_real_T(sp, c_u, i1);
  loop_ub = u->relativePositions->size[0] * u->relativePositions->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    c_u->data[i1] = u->relativePositions->data[i1];
  }

  b_y = NULL;
  m0 = emlrtCreateNumericArray(2, *(int32_T (*)[2])c_u->size, mxDOUBLE_CLASS,
    mxREAL);
  b_pData = emlrtMxGetPr(m0);
  i1 = 0;
  emxFree_real_T(sp, &c_u);
  for (loop_ub = 0; loop_ub < 3; loop_ub++) {
    for (i = 0; i < u->relativePositions->size[0]; i++) {
      b_pData[i1] = u->relativePositions->data[i + u->relativePositions->size[0]
        * loop_ub];
      i1++;
    }
  }

  emxInit_boolean_T(sp, &d_u, 2, true);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "relativePositions", b_y, 11);
  i1 = d_u->size[0] * d_u->size[1];
  d_u->size[0] = u->objectStatus->size[0];
  d_u->size[1] = 4;
  emxEnsureCapacity_boolean_T(sp, d_u, i1);
  loop_ub = u->objectStatus->size[0] * u->objectStatus->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    d_u->data[i1] = u->objectStatus->data[i1];
  }

  b_y = NULL;
  m0 = emlrtCreateLogicalArray(2, *(int32_T (*)[2])d_u->size);
  emlrtInitLogicalArray(u->objectStatus->size[0] << 2, m0, d_u->data);
  emlrtAssign(&b_y, m0);
  emlrtSetFieldR2017b(y, 0, "objectStatus", b_y, 12);
  emxFree_boolean_T(sp, &d_u);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_boolean_T *emxArray
 *                int32_T oldNumel
 * Return Type  : void
 */
static void emxEnsureCapacity_boolean_T(const emlrtStack *sp, emxArray_boolean_T
  *emxArray, int32_T oldNumel)
{
  int32_T newNumel;
  int32_T i;
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

    newData = emlrtCallocMex((uint32_T)i, sizeof(boolean_T));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(boolean_T) * oldNumel);
      if (emxArray->canFreeData) {
        emlrtFreeMex2018a(sp, (void *)emxArray->data);
      }
    }

    emxArray->data = (boolean_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_char_T *emxArray
 *                int32_T oldNumel
 * Return Type  : void
 */
static void emxEnsureCapacity_char_T(const emlrtStack *sp, emxArray_char_T
  *emxArray, int32_T oldNumel)
{
  int32_T newNumel;
  int32_T i;
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

    newData = emlrtCallocMex((uint32_T)i, sizeof(char_T));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(char_T) * oldNumel);
      if (emxArray->canFreeData) {
        emlrtFreeMex2018a(sp, (void *)emxArray->data);
      }
    }

    emxArray->data = (char_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_real_T *emxArray
 *                int32_T oldNumel
 * Return Type  : void
 */
static void emxEnsureCapacity_real_T(const emlrtStack *sp, emxArray_real_T
  *emxArray, int32_T oldNumel)
{
  int32_T newNumel;
  int32_T i;
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

    newData = emlrtCallocMex((uint32_T)i, sizeof(real_T));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(real_T) * oldNumel);
      if (emxArray->canFreeData) {
        emlrtFreeMex2018a(sp, (void *)emxArray->data);
      }
    }

    emxArray->data = (real_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_struct1_T *emxArray
 *                int32_T oldNumel
 * Return Type  : void
 */
static void emxEnsureCapacity_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  *emxArray, int32_T oldNumel)
{
  int32_T newNumel;
  int32_T i;
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

    newData = emlrtCallocMex((uint32_T)i, sizeof(struct1_T));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(struct1_T) * oldNumel);
      if (emxArray->canFreeData) {
        emlrtFreeMex2018a(sp, (void *)emxArray->data);
      }
    }

    emxArray->data = (struct1_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }

  if (oldNumel > newNumel) {
    emxTrim_struct1_T(sp, emxArray, newNumel, oldNumel);
  } else {
    if (oldNumel < newNumel) {
      emxExpand_struct1_T(sp, emxArray, oldNumel, newNumel);
    }
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_uint8_T *emxArray
 *                int32_T oldNumel
 * Return Type  : void
 */
static void emxEnsureCapacity_uint8_T(const emlrtStack *sp, emxArray_uint8_T
  *emxArray, int32_T oldNumel)
{
  int32_T newNumel;
  int32_T i;
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

    newData = emlrtCallocMex((uint32_T)i, sizeof(uint8_T));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(uint8_T) * oldNumel);
      if (emxArray->canFreeData) {
        emlrtFreeMex2018a(sp, (void *)emxArray->data);
      }
    }

    emxArray->data = (uint8_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_struct1_T *emxArray
 *                int32_T fromIndex
 *                int32_T toIndex
 * Return Type  : void
 */
static void emxExpand_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  *emxArray, int32_T fromIndex, int32_T toIndex)
{
  int32_T i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_struct1_T(sp, &emxArray->data[i], false);
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                struct0_T *pStruct
 * Return Type  : void
 */
static void emxFreeStruct_struct0_T(const emlrtStack *sp, struct0_T *pStruct)
{
  emxFree_struct1_T(sp, &pStruct->OBJECTS);
  emxFreeStruct_struct2_T(sp, &pStruct->TIME);
  emxFree_uint8_T(sp, &pStruct->globalIDvector);
  emxFree_char_T(sp, &pStruct->outputPath);
  emxFree_char_T(sp, &pStruct->phase);
  emxFree_char_T(sp, &pStruct->systemFile);
  emxFree_char_T(sp, &pStruct->outputFile);
}

/*
 * Arguments    : const emlrtStack *sp
 *                struct1_T *pStruct
 * Return Type  : void
 */
static void emxFreeStruct_struct1_T(const emlrtStack *sp, struct1_T *pStruct)
{
  emxFree_char_T(sp, &pStruct->name);
  emxFree_char_T(sp, &pStruct->class);
  emxFree_char_T(sp, &pStruct->symbol);
  emxFree_real_T(sp, &pStruct->relativePositions);
  emxFree_boolean_T(sp, &pStruct->objectStatus);
}

/*
 * Arguments    : const emlrtStack *sp
 *                struct2_T *pStruct
 * Return Type  : void
 */
static void emxFreeStruct_struct2_T(const emlrtStack *sp, struct2_T *pStruct)
{
  emxFree_real_T(sp, &pStruct->timeVector);
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_boolean_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_boolean_T(const emlrtStack *sp, emxArray_boolean_T
  **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if (((*pEmxArray)->data != (boolean_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->data);
    }

    emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->size);
    emlrtFreeMex2018a(sp, (void *)*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_char_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_char_T(const emlrtStack *sp, emxArray_char_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_char_T *)NULL) {
    if (((*pEmxArray)->data != (char_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->data);
    }

    emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->size);
    emlrtFreeMex2018a(sp, (void *)*pEmxArray);
    *pEmxArray = (emxArray_char_T *)NULL;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_real_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (real_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->data);
    }

    emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->size);
    emlrtFreeMex2018a(sp, (void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_struct1_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  **pEmxArray)
{
  int32_T numEl;
  int32_T i;
  if (*pEmxArray != (emxArray_struct1_T *)NULL) {
    if ((*pEmxArray)->data != (struct1_T *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }

      for (i = 0; i < numEl; i++) {
        emxFreeStruct_struct1_T(sp, &(*pEmxArray)->data[i]);
      }

      if ((*pEmxArray)->canFreeData) {
        emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->data);
      }
    }

    emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->size);
    emlrtFreeMex2018a(sp, (void *)*pEmxArray);
    *pEmxArray = (emxArray_struct1_T *)NULL;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_uint8_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_uint8_T(const emlrtStack *sp, emxArray_uint8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_uint8_T *)NULL) {
    if (((*pEmxArray)->data != (uint8_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->data);
    }

    emlrtFreeMex2018a(sp, (void *)(*pEmxArray)->size);
    emlrtFreeMex2018a(sp, (void *)*pEmxArray);
    *pEmxArray = (emxArray_uint8_T *)NULL;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                struct0_T *pStruct
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInitStruct_struct0_T(const emlrtStack *sp, struct0_T *pStruct,
  boolean_T doPush)
{
  emxInit_struct1_T(sp, &pStruct->OBJECTS, 2, doPush);
  emxInitStruct_struct2_T(sp, &pStruct->TIME, doPush);
  emxInit_uint8_T(sp, &pStruct->globalIDvector, 2, doPush);
  emxInit_char_T(sp, &pStruct->outputPath, 2, doPush);
  emxInit_char_T(sp, &pStruct->phase, 2, doPush);
  emxInit_char_T(sp, &pStruct->systemFile, 2, doPush);
  emxInit_char_T(sp, &pStruct->outputFile, 2, doPush);
}

/*
 * Arguments    : const emlrtStack *sp
 *                struct1_T *pStruct
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInitStruct_struct1_T(const emlrtStack *sp, struct1_T *pStruct,
  boolean_T doPush)
{
  emxInit_char_T(sp, &pStruct->name, 2, doPush);
  emxInit_char_T(sp, &pStruct->class, 2, doPush);
  emxInit_char_T(sp, &pStruct->symbol, 2, doPush);
  emxInit_real_T(sp, &pStruct->relativePositions, 2, doPush);
  emxInit_boolean_T(sp, &pStruct->objectStatus, 2, doPush);
}

/*
 * Arguments    : const emlrtStack *sp
 *                struct2_T *pStruct
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInitStruct_struct2_T(const emlrtStack *sp, struct2_T *pStruct,
  boolean_T doPush)
{
  emxInit_real_T(sp, &pStruct->timeVector, 2, doPush);
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_boolean_T **pEmxArray
 *                int32_T numDimensions
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInit_boolean_T(const emlrtStack *sp, emxArray_boolean_T
  **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
  emxArray_boolean_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_boolean_T *)emlrtMallocMex(sizeof(emxArray_boolean_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2018a(sp, (void *)pEmxArray, (void (*)(const
      void *, void *))emxFree_boolean_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (boolean_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_char_T **pEmxArray
 *                int32_T numDimensions
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInit_char_T(const emlrtStack *sp, emxArray_char_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_char_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_char_T *)emlrtMallocMex(sizeof(emxArray_char_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2018a(sp, (void *)pEmxArray, (void (*)(const
      void *, void *))emxFree_char_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (char_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_real_T **pEmxArray
 *                int32_T numDimensions
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)emlrtMallocMex(sizeof(emxArray_real_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2018a(sp, (void *)pEmxArray, (void (*)(const
      void *, void *))emxFree_real_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_struct1_T **pEmxArray
 *                int32_T numDimensions
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInit_struct1_T(const emlrtStack *sp, emxArray_struct1_T
  **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
  emxArray_struct1_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_struct1_T *)emlrtMallocMex(sizeof(emxArray_struct1_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2018a(sp, (void *)pEmxArray, (void (*)(const
      void *, void *))emxFree_struct1_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (struct1_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_uint8_T **pEmxArray
 *                int32_T numDimensions
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInit_uint8_T(const emlrtStack *sp, emxArray_uint8_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_uint8_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_uint8_T *)emlrtMallocMex(sizeof(emxArray_uint8_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2018a(sp, (void *)pEmxArray, (void (*)(const
      void *, void *))emxFree_uint8_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (uint8_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_struct1_T *emxArray
 *                int32_T fromIndex
 *                int32_T toIndex
 * Return Type  : void
 */
static void emxTrim_struct1_T(const emlrtStack *sp, emxArray_struct1_T *emxArray,
  int32_T fromIndex, int32_T toIndex)
{
  int32_T i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_struct1_T(sp, &emxArray->data[i]);
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real32_T y[3]
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[3])
{
  x_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                emxArray_real_T *ret
 * Return Type  : void
 */
static void fb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { 1, -1 };

  const boolean_T bv4[2] = { false, true };

  int32_T iv6[2];
  int32_T i9;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv4[0],
    iv6);
  i9 = ret->size[0] * ret->size[1];
  ret->size[0] = iv6[0];
  ret->size[1] = iv6[1];
  emxEnsureCapacity_real_T(sp, ret, i9);
  emlrtImportArrayR2015b(sp, src, ret->data, 8, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = y_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                emxArray_uint8_T *ret
 * Return Type  : void
 */
static void gb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_uint8_T *ret)
{
  static const int32_T dims[2] = { 1, -1 };

  const boolean_T bv5[2] = { false, true };

  int32_T iv7[2];
  int32_T i10;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "uint8", false, 2U, dims, &bv5[0],
    iv7);
  i10 = ret->size[0] * ret->size[1];
  ret->size[0] = iv7[0];
  ret->size[1] = iv7[1];
  emxEnsureCapacity_uint8_T(sp, ret, i10);
  emlrtImportArrayR2015b(sp, src, ret->data, 1, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : boolean_T
 */
static boolean_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId)
{
  boolean_T y;
  y = ab_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[3]
 */
static real_T (*hb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[3]
{
  real_T (*ret)[3];
  static const int32_T dims[1] = { 3 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[3])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[10]
 * Return Type  : void
 */
  static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[10])
{
  bb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[4]
 */
static real_T (*ib_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4]
{
  real_T (*ret)[4];
  static const int32_T dims[1] = { 4 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[4])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[9]
 * Return Type  : void
 */
  static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[9])
{
  cb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  db_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                emxArray_boolean_T *y
 * Return Type  : void
 */
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_boolean_T *y)
{
  eb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                struct2_T *y
 * Return Type  : void
 */
static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct2_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[10] = { "duration", "dt", "startTime",
    "idleTimeOut", "endTime", "timeVector", "currentTime", "frequency",
    "numSteps", "currentStep" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 10, fieldNames, 0U, &dims);
  thisId.fIdentifier = "duration";
  y->duration = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    0, "duration")), &thisId);
  thisId.fIdentifier = "dt";
  y->dt = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 1,
    "dt")), &thisId);
  thisId.fIdentifier = "startTime";
  y->startTime = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    2, "startTime")), &thisId);
  thisId.fIdentifier = "idleTimeOut";
  y->idleTimeOut = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 3, "idleTimeOut")), &thisId);
  thisId.fIdentifier = "endTime";
  y->endTime = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4,
    "endTime")), &thisId);
  thisId.fIdentifier = "timeVector";
  n_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 5,
    "timeVector")), &thisId, y->timeVector);
  thisId.fIdentifier = "currentTime";
  y->currentTime = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 6, "currentTime")), &thisId);
  thisId.fIdentifier = "frequency";
  y->frequency = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    7, "frequency")), &thisId);
  thisId.fIdentifier = "numSteps";
  y->numSteps = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    8, "numSteps")), &thisId);
  thisId.fIdentifier = "currentStep";
  y->currentStep = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 9, "currentStep")), &thisId);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  fb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                emxArray_uint8_T *y
 * Return Type  : void
 */
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_uint8_T *y)
{
  gb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *objectID
 *                const char_T *identifier
 * Return Type  : uint8_T
 */
static uint8_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *objectID,
  const char_T *identifier)
{
  uint8_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(objectID), &thisId);
  emlrtDestroyArray(&objectID);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *globalVelocity_k
 *                const char_T *identifier
 * Return Type  : real_T (*)[3]
 */
static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *globalVelocity_k, const char_T *identifier))[3]
{
  real_T (*y)[3];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = r_emlrt_marshallIn(sp, emlrtAlias(globalVelocity_k), &thisId);
  emlrtDestroyArray(&globalVelocity_k);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[3]
 */
  static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[3]
{
  real_T (*y)[3];
  y = hb_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *quaternion_k
 *                const char_T *identifier
 * Return Type  : real_T (*)[4]
 */
static real_T (*s_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *quaternion_k, const char_T *identifier))[4]
{
  real_T (*y)[4];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = t_emlrt_marshallIn(sp, emlrtAlias(quaternion_k), &thisId);
  emlrtDestroyArray(&quaternion_k);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[4]
 */
  static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[4]
{
  real_T (*y)[4];
  y = ib_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *idleStatus_k
 *                const char_T *identifier
 * Return Type  : boolean_T
 */
static boolean_T u_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *idleStatus_k, const char_T *identifier)
{
  boolean_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(idleStatus_k), &thisId);
  emlrtDestroyArray(&idleStatus_k);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : uint8_T
 */
static uint8_T v_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  uint8_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "uint8", false, 0U, &dims);
  ret = *(uint8_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                emxArray_char_T *ret
 * Return Type  : void
 */
static void w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_char_T *ret)
{
  static const int32_T dims[2] = { 1, -1 };

  const boolean_T bv1[2] = { false, true };

  int32_T iv3[2];
  int32_T i2;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims, &bv1[0],
    iv3);
  i2 = ret->size[0] * ret->size[1];
  ret->size[0] = iv3[0];
  ret->size[1] = iv3[1];
  emxEnsureCapacity_char_T(sp, ret, i2);
  emlrtImportArrayR2015b(sp, src, ret->data, 1, false);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real32_T ret[3]
 * Return Type  : void
 */
static void x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[3])
{
  static const int32_T dims[2] = { 1, 3 };

  int32_T i3;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 2U, dims);
  for (i3 = 0; i3 < 3; i3++) {
    ret[i3] = (*(real32_T (*)[3])emlrtMxGetData(src))[i3];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const mxArray * const prhs[5]
 *                int32_T nlhs
 *                const mxArray *plhs[1]
 * Return Type  : void
 */
void OMAS_updateGlobalStates_api(const mxArray * const prhs[5], int32_T nlhs,
  const mxArray *plhs[1])
{
  struct0_T SIM;
  struct1_T METAObjUpdate;
  uint8_T objectID;
  real_T (*globalVelocity_k)[3];
  real_T (*quaternion_k)[4];
  boolean_T idleStatus_k;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInitStruct_struct0_T(&st, &SIM, true);
  emxInitStruct_struct1_T(&st, &METAObjUpdate, true);

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "SIM", &SIM);
  objectID = p_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "objectID");
  globalVelocity_k = q_emlrt_marshallIn(&st, emlrtAlias(prhs[2]),
    "globalVelocity_k");
  quaternion_k = s_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "quaternion_k");
  idleStatus_k = u_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "idleStatus_k");

  /* Invoke the target function */
  OMAS_updateGlobalStates(&SIM, objectID, *globalVelocity_k, *quaternion_k,
    idleStatus_k, &METAObjUpdate);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(&st, &METAObjUpdate);
  emxFreeStruct_struct1_T(&st, &METAObjUpdate);
  emxFreeStruct_struct0_T(&st, &SIM);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void OMAS_updateGlobalStates_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  OMAS_updateGlobalStates_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void OMAS_updateGlobalStates_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void OMAS_updateGlobalStates_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_OMAS_updateGlobalStates_api.c
 *
 * [EOF]
 */
