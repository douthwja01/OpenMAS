/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Jan-2019 16:40:18
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
#include "OMAS_updateGlobalStates.h"
#include "main.h"
#include "OMAS_updateGlobalStates_terminate.h"
#include "OMAS_updateGlobalStates_emxAPI.h"
#include "OMAS_updateGlobalStates_initialize.h"

/* Function Declarations */
static void argInit_10x1_real_T(double result[10]);
static void argInit_1x3_real32_T(float result[3]);
static emxArray_char_T *argInit_1xUnbounded_char_T(void);
static emxArray_real_T *argInit_1xUnbounded_real_T(void);
static emxArray_struct1_T *argInit_1xUnbounded_struct1_T(void);
static emxArray_uint8_T *argInit_1xUnbounded_uint8_T(void);
static void argInit_3x1_real_T(double result[3]);
static void argInit_3x3_real_T(double result[9]);
static void argInit_4x1_real_T(double result[4]);
static emxArray_real_T *argInit_Unboundedx3_real_T(void);
static emxArray_boolean_T *argInit_Unboundedx4_boolean_T(void);
static boolean_T argInit_boolean_T(void);
static char argInit_char_T(void);
static float argInit_real32_T(void);
static double argInit_real_T(void);
static void argInit_struct0_T(struct0_T *result);
static void argInit_struct1_T(struct1_T *result);
static void argInit_struct2_T(struct2_T *result);
static unsigned char argInit_uint8_T(void);
static void main_OMAS_updateGlobalStates(void);

/* Function Definitions */

/*
 * Arguments    : double result[10]
 * Return Type  : void
 */
static void argInit_10x1_real_T(double result[10])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 10; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : float result[3]
 * Return Type  : void
 */
static void argInit_1x3_real32_T(float result[3])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 3; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real32_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : emxArray_char_T *
 */
static emxArray_char_T *argInit_1xUnbounded_char_T(void)
{
  emxArray_char_T *result;
  static int iv1[2] = { 1, 2 };

  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_char_T(2, iv1);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[result->size[0] * idx1] = argInit_char_T();
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : emxArray_real_T *
 */
static emxArray_real_T *argInit_1xUnbounded_real_T(void)
{
  emxArray_real_T *result;
  static int iv4[2] = { 1, 2 };

  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(2, iv4);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[result->size[0] * idx1] = argInit_real_T();
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : emxArray_struct1_T *
 */
static emxArray_struct1_T *argInit_1xUnbounded_struct1_T(void)
{
  emxArray_struct1_T *result;
  static int iv0[2] = { 1, 2 };

  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_struct1_T(2, iv0);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    argInit_struct1_T(&result->data[result->size[0] * idx1]);
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : emxArray_uint8_T *
 */
static emxArray_uint8_T *argInit_1xUnbounded_uint8_T(void)
{
  emxArray_uint8_T *result;
  static int iv5[2] = { 1, 2 };

  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_uint8_T(2, iv5);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[result->size[0] * idx1] = argInit_uint8_T();
  }

  return result;
}

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[9]
 * Return Type  : void
 */
static void argInit_3x3_real_T(double result[9])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : double result[4]
 * Return Type  : void
 */
static void argInit_4x1_real_T(double result[4])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : emxArray_real_T *
 */
static emxArray_real_T *argInit_Unboundedx3_real_T(void)
{
  emxArray_real_T *result;
  static int iv2[2] = { 2, 3 };

  int idx0;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(2, iv2);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : emxArray_boolean_T *
 */
static emxArray_boolean_T *argInit_Unboundedx4_boolean_T(void)
{
  emxArray_boolean_T *result;
  static int iv3[2] = { 2, 4 };

  int idx0;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_boolean_T(2, iv3);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < 4; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_boolean_T();
    }
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : boolean_T
 */
static boolean_T argInit_boolean_T(void)
{
  return false;
}

/*
 * Arguments    : void
 * Return Type  : char
 */
static char argInit_char_T(void)
{
  return '?';
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
 * Arguments    : struct0_T *result
 * Return Type  : void
 */
static void argInit_struct0_T(struct0_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->OBJECTS = argInit_1xUnbounded_struct1_T();
  argInit_struct2_T(&result->TIME);
  result->conditionTolerance = argInit_real_T();
  result->globalIDvector = argInit_1xUnbounded_uint8_T();
  result->gui = argInit_boolean_T();
  result->monteCarloMode = argInit_boolean_T();
  result->outputPath = argInit_1xUnbounded_char_T();
  result->phase = argInit_1xUnbounded_char_T();
  result->systemFile = argInit_1xUnbounded_char_T();
  result->threadPool = argInit_boolean_T();
  result->totalAgents = argInit_uint8_T();
  result->totalMiscs = argInit_uint8_T();
  result->totalObjects = argInit_uint8_T();
  result->totalObstacles = argInit_uint8_T();
  result->totalWaypoints = argInit_uint8_T();
  result->verbosity = argInit_uint8_T();
  result->visabilityModifier = argInit_real_T();
  result->warningDistance = argInit_real_T();
  result->outputFile = argInit_1xUnbounded_char_T();
}

/*
 * Arguments    : struct1_T *result
 * Return Type  : void
 */
static void argInit_struct1_T(struct1_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->objectID = argInit_uint8_T();
  result->name = argInit_1xUnbounded_char_T();
  result->class = argInit_1xUnbounded_char_T();
  result->type = argInit_uint8_T();
  argInit_1x3_real32_T(result->colour);
  result->symbol = argInit_1xUnbounded_char_T();
  result->radius = argInit_real_T();
  result->detectionRadius = argInit_real_T();
  result->idleStatus = argInit_boolean_T();
  argInit_10x1_real_T(result->globalState);
  argInit_3x3_real_T(result->R);
  result->relativePositions = argInit_Unboundedx3_real_T();
  result->objectStatus = argInit_Unboundedx4_boolean_T();
}

/*
 * Arguments    : struct2_T *result
 * Return Type  : void
 */
static void argInit_struct2_T(struct2_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->duration = argInit_real_T();
  result->dt = argInit_real_T();
  result->startTime = argInit_real_T();
  result->idleTimeOut = argInit_real_T();
  result->endTime = argInit_real_T();
  result->timeVector = argInit_1xUnbounded_real_T();
  result->currentTime = argInit_real_T();
  result->frequency = argInit_real_T();
  result->numSteps = argInit_real_T();
  result->currentStep = argInit_real_T();
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
 * Return Type  : void
 */
static void main_OMAS_updateGlobalStates(void)
{
  struct1_T METAObjUpdate;
  struct0_T SIM;
  double dv0[3];
  double dv1[4];
  emxInit_struct1_T(&METAObjUpdate);

  /* Initialize function 'OMAS_updateGlobalStates' input arguments. */
  /* Initialize function input argument 'SIM'. */
  argInit_struct0_T(&SIM);

  /* Initialize function input argument 'globalVelocity_k'. */
  /* Initialize function input argument 'quaternion_k'. */
  /* Call the entry-point 'OMAS_updateGlobalStates'. */
  argInit_3x1_real_T(dv0);
  argInit_4x1_real_T(dv1);
  OMAS_updateGlobalStates(&SIM, argInit_uint8_T(), dv0, dv1, argInit_boolean_T(),
    &METAObjUpdate);
  emxDestroy_struct1_T(METAObjUpdate);
  emxDestroy_struct0_T(SIM);
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
  OMAS_updateGlobalStates_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_OMAS_updateGlobalStates();

  /* Terminate the application.
     You do not need to do this more than one time. */
  OMAS_updateGlobalStates_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
