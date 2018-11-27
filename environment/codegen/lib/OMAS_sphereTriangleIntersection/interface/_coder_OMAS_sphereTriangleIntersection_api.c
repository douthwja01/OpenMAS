/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_OMAS_sphereTriangleIntersection_api.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 07-Nov-2018 17:03:23
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_OMAS_sphereTriangleIntersection_api.h"
#include "_coder_OMAS_sphereTriangleIntersection_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131466U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "OMAS_sphereTriangleIntersection",   /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[3];
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *r, const
  char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[3];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *P, const
  char_T *identifier))[3];
static const mxArray *emlrt_marshallOut(const boolean_T u);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[3]
 */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[3]
{
  real_T (*y)[3];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *r
 *                const char_T *identifier
 * Return Type  : real_T
 */
  static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *r, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(r), &thisId);
  emlrtDestroyArray(&r);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[3]
 */
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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
 *                const mxArray *P
 *                const char_T *identifier
 * Return Type  : real_T (*)[3]
 */
  static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *P, const
  char_T *identifier))[3]
{
  real_T (*y)[3];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(P), &thisId);
  emlrtDestroyArray(&P);
  return y;
}

/*
 * Arguments    : const boolean_T u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const boolean_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateLogicalScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
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
void OMAS_sphereTriangleIntersection_api(const mxArray * const prhs[5], int32_T
  nlhs, const mxArray *plhs[1])
{
  const mxArray *prhs_copy_idx_2;
  const mxArray *prhs_copy_idx_3;
  const mxArray *prhs_copy_idx_4;
  real_T (*P)[3];
  real_T r;
  real_T (*A)[3];
  real_T (*B)[3];
  real_T (*C)[3];
  boolean_T haveCollided;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, false, -1);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, false, -1);

  /* Marshall function inputs */
  P = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "P");
  r = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "r");
  A = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "A");
  B = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_3), "B");
  C = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_4), "C");

  /* Invoke the target function */
  haveCollided = OMAS_sphereTriangleIntersection(*P, r, *A, *B, *C);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(haveCollided);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void OMAS_sphereTriangleIntersection_atexit(void)
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
  OMAS_sphereTriangleIntersection_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void OMAS_sphereTriangleIntersection_initialize(void)
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
void OMAS_sphereTriangleIntersection_terminate(void)
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
 * File trailer for _coder_OMAS_sphereTriangleIntersection_api.c
 *
 * [EOF]
 */
