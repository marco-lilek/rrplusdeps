
/*

      joinvect.h : Function prototypes and data types for the standalone
               "joinvect" application. 

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

#ifdef HAVE_STDBOOL_H
#include<stdbool.h>
#endif

#ifdef HAVE_STRING_H
#include<string.h>
#endif

#include"nplCurve.h"
#include"argtable2.h"
#include"../utilib/mangle.h"
#include"../utilib/ordie.h"

/* progressbar.c */

void init_progressbar(int nevents);
void update_progressbar();

/* tsporder.c */

void tsp_recurser(int level,int ncurves,bool *used,double thislen,double *minlen,
		  int order[],int workingorder[]);
void tsp_order(int ncurves,nplCurve **curves,int *order,double twidth);
void greedy_order(int ncurves,nplCurve **curves,int *order,double twidth);

nplc_vector nplc_start(nplCurve *L,int cp);
nplc_vector nplc_end(nplCurve *L,int cp);
nplc_vector nplc_sf(nplCurve *L,int cp,int sf);

/* torusdistance.c */

/* Finds the longest edge in an nplCurve */
void longest_edge(const nplCurve * const L, int *cp, int *vt0, int *vt1, double *len,
		  double twidth);
double torus_distance(nplc_vector A, nplc_vector B,double twidth);
