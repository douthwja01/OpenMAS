/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_geometry.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "OMAS_updateGlobalStates.h"
#include "OMAS_geometry.h"

/* Function Definitions */

/*
 * Arguments    : double q[4]
 *                double R[9]
 * Return Type  : void
 */
void c_OMAS_geometry_quaternionToRot(double q[4], double R[9])
{
  int k;
  double y;
  double z1[4];
  for (k = 0; k < 4; k++) {
    z1[k] = q[k] * q[k];
  }

  y = z1[0];
  for (k = 0; k < 3; k++) {
    y += z1[k + 1];
  }

  y = sqrt(y);
  for (k = 0; k < 4; k++) {
    q[k] /= y;
  }

  R[0] = ((q[0] * q[0] + q[1] * q[1]) - q[2] * q[2]) - q[3] * q[3];
  R[3] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  R[6] = 2.0 * (q[0] * q[2] + q[1] * q[3]);
  R[1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  R[4] = ((q[0] * q[0] - q[1] * q[1]) + q[2] * q[2]) - q[3] * q[3];
  R[7] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
  R[2] = 2.0 * (q[1] * q[3] - q[0] * q[2]);
  R[5] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
  R[8] = ((q[0] * q[0] - q[1] * q[1]) - q[2] * q[2]) + q[3] * q[3];
}

/*
 * File trailer for OMAS_geometry.c
 *
 * [EOF]
 */
