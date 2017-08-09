/*
 *  cdoublemodel.c

 This file, part of writhecomp, includes standard C doubles as an arithmetic
 model for writhe computations. We will also probably include a "long double"
 model, to see what effect this has on the computations. 

 */

#include "writhecomp.h"

void init_cdouble_model(int precision) 

     /* Procedure just initializes the standard constants np_PI, np_MINUS. np_EPS. */

{

  /* First, we need memory. */

  np_PI = (NUMPTR) calloc(1,sizeof(double));
  np_MINUS = (NUMPTR) calloc(1,sizeof(double));
  np_EPS = (NUMPTR) calloc(1,sizeof(double));
  np_ZERO = (NUMPTR) calloc(1,sizeof(double));

  /* Now, we add values. */

  *(double *)(np_PI) = 3.1415926535897932;
  *(double *)(np_MINUS) = -1.0;
  *(double *)(np_EPS) = 1e-11;
  *(double *)(np_ZERO) = 0.0;

}

void kill_cdouble_model()

     /* Procedure eliminates the standard constants. */

{
  free(np_PI);
  free(np_MINUS);
  free(np_EPS);
  free(np_ZERO);
}

NUMPTR double_to_cdouble(double x)

     /* Procedure converts a double to a NUMPTR. */

{
  NUMPTR ret;

  ret = calloc(1,sizeof(double));
  *((double *)(ret)) = x;

  return ret;
}

double cdouble_to_double(NUMPTR x)

     /* Procedure converts a NUMPTR to a cdouble. */

{
  return *((double *)(x));
}

char *cdouble_to_string(NUMPTR x)

     /* Procedure converts a double to a string. */

{
  char *result;

  result = calloc(25,sizeof(char));
  sprintf(result,"%.15e",*(double *)(x));

  return result;
}

NUMPTR cdouble_add(NUMPTR a,NUMPTR b)

     /* Procedure adds doubles, with some annoying memory overhead. */

{
  NUMPTR result;

  result = (NUMPTR)(calloc(1,sizeof(double)));
  *((double *)(result)) = *(double *)(a) + *(double *)(b);
  
  return result;
}

NUMPTR cdouble_mul(NUMPTR a,NUMPTR b)

     /* Procedure multiplies doubles, with some annoying overhead. */

{
  NUMPTR result;

  result = (NUMPTR)(calloc(1,sizeof(double)));
  *((double *)(result)) = (*(double *)(a)) * (*(double *)(b));
  
  return result;
}

NUMPTR cdouble_div(NUMPTR a,NUMPTR b)

     /* Procedure divides doubles, with overhead. */

{
  NUMPTR result;

  result = (NUMPTR)(calloc(1,sizeof(double)));
  *((double *)(result)) = *(double *)(a) / *(double *)(b);
  
  return result;
}

NUMPTR cdouble_acos(NUMPTR theta)

     /* Procedure takes the arc cosine of theta using the standard 
	library function. Note that the numerical quality of this 
	result may vary greatly from machine to machine, so this is
	not considered a stable operation. */

     /* A future release of writhecomp will probably replace this 
	with better library code from some source. */

{
  NUMPTR result;

  result = (NUMPTR)(calloc(1,sizeof(double)));

  if (*(double *)(theta) > 1.0) {

    *(double *)(result) = 0.0;

  } else if (*(double *)(theta) < -1.0) {

    *(double *)(result) = 3.141592653589793;

  } else {

    *((double *)(result)) = acos(*(double *)(theta));

  }

  return result;
}

NUMPTR cdouble_sqrt(NUMPTR x)

     /* Procedure takes the square root of x, using the standard library function. */

{
  NUMPTR result;

  result = (NUMPTR)(calloc(1,sizeof(double)));
  *((double *)(result)) = sqrt(*(double *)(x));
  
  return result;
}

int cdouble_leq(NUMPTR a, NUMPTR b) 

     /* Procedure compares doubles. */

{
  if (*(double *)(a) <= *(double *)(b)) {

    return TRUE;

  } else {

    return FALSE;

  }
}

