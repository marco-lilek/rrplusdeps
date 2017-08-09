/*************************************************************/  
/*                                                           */
/*  writhecomp.h : Master header file for the writhecomp     */
/*                 project. Note that all new arithmetic     */
/*                 models must declare their stub functions  */
/*                 below, and include this header.           */
/*                                                           */
/*************************************************************/

#ifndef __WRITHECOMP_H__
#define __WRITHECOMP_H__

#if     HAVE_CONFIG_H
#include<config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#ifdef HAVE_MATH_H
#include <math.h>
#endif

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif

#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#include <argtable2.h>

#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif

#ifdef PARI
#include "genpari.h"
#endif

#include "m_apm.h"  // This is a local file which is now distributed with writhecomp.

#include "plCurve.h" // We depend on the plCurve library.

//#include "../libraries/linearalgebra.h"
//#include "../libraries/polylines.h"

#define DBL_PRECISION 10        
/* The number of decimal digits of precision in a C double. */

double  timing;			
/* The number of ms required to compute one pair_writhe + crossing_number. */

FILE	*log_file;

extern int TIMING;         /* Just estimate the time for a computation? no. */
extern int SAFE;	   /* Safe mode self-checking?  Default- no.*/ 
extern int VERBOSE;	   /* Report everything to the user? Default- no. */
extern int DEC_PRECISION;  /* Digits of (decimal) precision to use. */

#ifdef PARI

int     QD_VERTS = {250};

/* Pari measures precision in long words, while we measure precision
   in decimal digits. Setting PRECISION to the Pari defined constant
   DEFAULTPREC guarantees 19 digits, no matter whether we are running
   on a 32 or 64-bit architecture. */

/* Later, we will set DEC_PRECISION ourselves on the command line, and
   convert to a Pari precision in the code. */  		     

#endif

unsigned long int PARISTACK;

#ifdef PARI
#include"parimodel.h"
#endif

#ifndef FALSE 
#define FALSE (1 == 0)
#endif

#ifndef TRUE
#define TRUE (1 == 1)
#endif

typedef void *NUMPTR;
typedef NUMPTR (NUMVECTOR)[3];

int arithmodel_ok();

/********  Function Prototypes for the MAPM model *********/

#include"m_apm.h"

extern int MAPM_running;       /* This global tells us whether MAPM is already
				  running. This is useful if we bring up and shut
				  down the the model repeatedly, as in SAFE mode. */

void   init_mapm_model(int precision);
void   kill_mapm_model();
void   kill_mapm(NUMPTR x);

NUMPTR double_to_mapm(double x);
double mapm_to_double(NUMPTR x);
char  *mapm_to_string(NUMPTR x);

NUMPTR mapm_add(NUMPTR a,NUMPTR b);
NUMPTR mapm_mul(NUMPTR a,NUMPTR b);
NUMPTR mapm_div(NUMPTR a,NUMPTR b);

NUMPTR mapm_acos(NUMPTR theta);
NUMPTR mapm_sqrt(NUMPTR x);
int    mapm_leq(NUMPTR a,NUMPTR b);

/*********  Function Prototypes for the C Double model *****/

void   init_cdouble_model(int precision);
void   kill_cdouble_model();

NUMPTR double_to_cdouble(double x);
double cdouble_to_double(NUMPTR x);
char  *cdouble_to_string(NUMPTR x);

NUMPTR cdouble_add(NUMPTR a,NUMPTR b);
NUMPTR cdouble_mul(NUMPTR a,NUMPTR b);
NUMPTR cdouble_div(NUMPTR a,NUMPTR b);

NUMPTR cdouble_acos(NUMPTR theta);
NUMPTR cdouble_sqrt(NUMPTR x);
int    cdouble_leq(NUMPTR a,NUMPTR b);

/******************************************************************************/
/* We now define a set of standard (core) arithmetic model functions which    */
/* we will use for the computation. Everything is built carefully in terms of */
/* these guys. The intent is to make it easy to add new arithmetic models in  */
/* the future while preserving Banchoff's algorithm in the code.              */
/******************************************************************************/

NUMPTR  np_PI;     /* We expect these constants to be set by init_arith_model. */
NUMPTR  np_MINUS;  /* And freed by kill_arith_model. */
NUMPTR  np_EPS;    
NUMPTR  np_ZERO;

void   (*init_arith_model)(int precision);
void   (*kill_arith_model)();
void   (*kill_np)(NUMPTR x);

NUMPTR (*double_to_np)(double x);
double (*np_to_double)(NUMPTR x);
char  *(*np_to_string)(NUMPTR x);

NUMPTR (*np_add)(NUMPTR a,NUMPTR b);
NUMPTR (*np_mul)(NUMPTR a,NUMPTR b);
NUMPTR (*np_div)(NUMPTR a,NUMPTR b);
NUMPTR (*np_acos)(NUMPTR theta);
NUMPTR (*np_sqrt)(NUMPTR x);
int    (*np_leq)(NUMPTR a, NUMPTR b);

/***************************************************************/
/* We now define the model-independant writhe functions.       */
/* These are stored in anymodel.c.                             */
/***************************************************************/

NUMPTR    np_plusequal(NUMPTR a,NUMPTR b);
NUMPTR    np_minusequal(NUMPTR a,NUMPTR b);
NUMPTR    *nv_diff(NUMVECTOR a,NUMVECTOR b);
NUMPTR    *nv_cross(NUMVECTOR a,NUMVECTOR b);
NUMPTR    nv_dot(NUMVECTOR a,NUMVECTOR b);
void      nv_normalize(NUMVECTOR a);
void      kill_nv(NUMVECTOR a);

NUMPTR    np_dihedral_sum(NUMVECTOR X_one,NUMVECTOR X_two,
                          NUMVECTOR Y_one,NUMVECTOR Y_two);

NUMPTR    np_writhe(plCurve *pline,int comp1,int comp2); 
// Computes pairwise write for components comp1 and comp2

double    writhe(plCurve *pline,int precision,double *pairvalues);
// Compute the writhe for pline in precision digits of precision. 
// If pairvalues is non-NULL, it is expected to be an array of length
// pline->nc * pline->nc. On return, it will contain the pairwise Gauss 
// integrals of the components of pline.

int       test_model();

/************************************************************/ 

double     ni_integrand[100][100];
double     np_integrand[100][100];

#endif
