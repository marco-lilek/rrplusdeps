/*

     orientvect.h : Function prototypes and data types for the standalone
      "orientvect" application. 

*/

#if HAVE_CONFIG_H
#include<config.h>
#endif

#ifdef HAVE_ASSERT_H
#include<assert.h>
#endif

#ifdef HAVE_MATH_H
#include<math.h>
#endif

#ifdef HAVE_STDLIB_H
#include<stdlib.h>
#endif

#ifdef HAVE_STRING_H
#include<string.h>
#endif

#ifdef HAVE_GSL_GSL_MATH_H
#include<gsl/gsl_math.h>
#endif

#ifdef HAVE_GSL_GSL_EIGEN_H
#include<gsl/gsl_eigen.h>
#endif

#include"plCurve.h"

extern int  octrope_error_num;
extern char octrope_error_str[80];


