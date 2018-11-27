/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Nov-2018 17:59:35
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "OMAS_getTrajectoryData.h"
#include "main.h"
#include "OMAS_getTrajectoryData_terminate.h"
#include "OMAS_getTrajectoryData_emxAPI.h"
#include "OMAS_getTrajectoryData_initialize.h"

/* Function Declarations */
static emxArray_uint8_T *argInit_1xUnbounded_uint8_T(void);
static float argInit_real32_T(void);
static double argInit_real_T(void);
static unsigned char argInit_uint8_T(void);
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void);
static void main_OMAS_getTrajectoryData(void);

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : emxArray_uint8_T *
 */
static emxArray_uint8_T *argInit_1xUnbounded_uint8_T(void)
{
  emxArray_uint8_T *result;
  static int iv1[2] = { 1, 2 };

  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_uint8_T(2, iv1);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[result->size[0] * idx1] = argInit_uint8_T();
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : float
 */
static float argInit_real32_T(void)
{
  return 0.0F;
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : unsigned char
 */
static unsigned char argInit_uint8_T(void)
{
  return 0U;
}

/*
 * Arguments    : void
 * Return Type  : emxArray_real_T *
 */
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void)
{
  emxArray_real_T *result;
  static int iv0[2] = { 2, 2 };

  int idx0;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(2, iv0);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_OMAS_getTrajectoryData(void)
{
  emxArray_real_T *objectStateData;
  emxArray_real_T *trajectoryMatrix;
  emxArray_uint8_T *globalIDvector;
  unsigned int indexOfLastState;
  emxInitArray_real_T(&objectStateData, 2);

  /* Initialize function 'OMAS_getTrajectoryData' input arguments. */
  /* Initialize function input argument 'trajectoryMatrix'. */
  trajectoryMatrix = c_argInit_UnboundedxUnbounded_r();

  /* Initialize function input argument 'globalIDvector'. */
  globalIDvector = argInit_1xUnbounded_uint8_T();

  /* Call the entry-point 'OMAS_getTrajectoryData'. */
  OMAS_getTrajectoryData(trajectoryMatrix, globalIDvector, argInit_uint8_T(),
    argInit_real32_T(), objectStateData, &indexOfLastState);
  emxDestroyArray_real_T(objectStateData);
  emxDestroyArray_uint8_T(globalIDvector);
  emxDestroyArray_real_T(trajectoryMatrix);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  OMAS_getTrajectoryData_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_OMAS_getTrajectoryData();

  /* Terminate the application.
     You do not need to do this more than one time. */
  OMAS_getTrajectoryData_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
