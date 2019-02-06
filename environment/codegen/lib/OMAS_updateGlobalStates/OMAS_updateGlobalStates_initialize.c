/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_updateGlobalStates_initialize.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "OMAS_updateGlobalStates.h"
#include "OMAS_updateGlobalStates_initialize.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void OMAS_updateGlobalStates_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * File trailer for OMAS_updateGlobalStates_initialize.c
 *
 * [EOF]
 */
