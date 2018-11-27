/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: OMAS_sphereTriangleIntersection.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 07-Nov-2018 17:03:23
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "OMAS_sphereTriangleIntersection.h"

/* Function Definitions */

/*
 * INPUT CHECK
 * Arguments    : const double P[3]
 *                double r
 *                double A[3]
 *                double B[3]
 *                double C[3]
 * Return Type  : boolean_T
 */
boolean_T OMAS_sphereTriangleIntersection(const double P[3], double r, double A
  [3], double B[3], double C[3])
{
  double rr;
  int i;
  double V[3];
  double AB[3];
  double BC[3];
  double d;
  double e;
  double aa;
  double ab;
  double ac;
  double bb;
  double bc;
  double cc;
  double d1;
  double d2;
  double d3;
  double e1;
  double e2;
  double e3;
  double c;
  double b_c;
  double c_c;
  double d_c;
  double e_c;
  double f_c;
  double b_AB;
  boolean_T isSeparated;
  double b_BC;
  double b_V;

  /*  This function is designed to determine whether a sphere collides with a */
  /*  give triangle */
  /*  CHECK PLANAR SEPARATION */
  rr = r * r;
  for (i = 0; i < 3; i++) {
    d = A[i] - P[i];
    e = B[i] - P[i];
    aa = C[i] - P[i];
    AB[i] = e - d;
    BC[i] = aa - d;
    A[i] = d;
    B[i] = e;
    C[i] = aa;
  }

  V[0] = AB[1] * BC[2] - AB[2] * BC[1];
  V[1] = AB[2] * BC[0] - AB[0] * BC[2];
  V[2] = AB[0] * BC[1] - AB[1] * BC[0];
  d = 0.0;
  e = 0.0;

  /*  CHECK VERTICES */
  aa = 0.0;
  ab = 0.0;
  ac = 0.0;
  bb = 0.0;
  bc = 0.0;
  cc = 0.0;

  /*  CHECK TRIANGLE EDGES */
  for (i = 0; i < 3; i++) {
    d += A[i] * V[i];
    e += V[i] * V[i];
    aa += A[i] * A[i];
    ab += A[i] * B[i];
    ac += A[i] * C[i];
    bb += B[i] * B[i];
    bc += B[i] * C[i];
    cc += C[i] * C[i];
    AB[i] = B[i] - A[i];
    BC[i] = C[i] - B[i];
    V[i] = A[i] - C[i];
  }

  d1 = ab - aa;
  d2 = bc - bb;
  d3 = ac - cc;
  e1 = 0.0;
  e2 = 0.0;
  e3 = 0.0;
  for (i = 0; i < 3; i++) {
    e1 += AB[i] * AB[i];
    e2 += BC[i] * BC[i];
    e3 += V[i] * V[i];
  }

  c = 0.0;
  b_c = 0.0;
  c_c = 0.0;
  d_c = 0.0;
  e_c = 0.0;
  f_c = 0.0;
  for (i = 0; i < 3; i++) {
    b_AB = A[i] * e1 - d1 * AB[i];
    b_BC = B[i] * e2 - d2 * BC[i];
    b_V = C[i] * e3 - d3 * V[i];
    c += b_AB * b_AB;
    b_c += b_AB * (C[i] * e1 - b_AB);
    c_c += b_BC * b_BC;
    d_c += b_BC * (A[i] * e2 - b_BC);
    e_c += b_V * b_V;
    f_c += b_V * (B[i] * e3 - b_V);
  }

  if ((d * d > rr * e) || ((aa > rr) && (ab > aa) && (ac > aa)) || ((bb > rr) &&
       (ab > bb) && (bc > bb)) || ((cc > rr) && (ac > cc) && (bc > cc)) || ((c >
        rr * e1 * e1) && (b_c > 0.0)) || ((c_c > rr * e2 * e2) && (d_c > 0.0)) ||
      ((e_c > rr * e3 * e3) && (f_c > 0.0))) {
    isSeparated = true;
  } else {
    isSeparated = false;
  }

  return !isSeparated;
}

/*
 * File trailer for OMAS_sphereTriangleIntersection.c
 *
 * [EOF]
 */
