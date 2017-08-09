/* ktcurve.c.  Generated by configure.  */
/*
 *
 * A library for generating plCurves from curvature and torsion using GSL.
 *
 *  $Id: ktcurve.c,v 1.1 2008/01/24 05:05:01 cantarel Exp $
 *
 */

/* Copyright 2008 The University of Georgia */

/* This file is part of vecttools.

vecttools is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

vecttools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vecttools; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

/* 
   Note: libktcurve requires the GSL library and the plCurve library.
*/


#include "ktcurve.h"
#include "plCurve.h"

#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#ifdef HAVE_MATH_H
  #include "math.h"
#endif

#ifdef HAVE_GSL_GSL_ERRNO_H
  #include "gsl/gsl_errno.h"
#endif

#ifdef HAVE_GSL_GSL_ODEIV_H
  #include "gsl/gsl_odeiv.h"
#endif

int  kterrcode;
char kterrstring[1024];

typedef struct {

  ktfunc *curvature;
  ktfunc *torsion;
  void   *curvedata;
  
} ktdata_type;

int kt_integration_func(double t, const double y[12], double dydt[12], void *params)

/* This function will be called by the GSL ODE solver. 
   It extracts the (12-dimensional!) system of ODEs which 
   define the Frenet-Serret equations from the initial values in y:
   
   gamma(t) = (y[0],y[1],y[2])
   T(t) = (y[3],y[4],y[5])
   N(t) = (y[6],y[7],y[8])
   B(t) = (y[9],y[10],y[11]). 
   
   and returns the values of 
   
   gamma'(t) = (dydt[0],dydt[1],dydt[2]),
   T'(t) = (dydt[3],dydt[4],dydt[5])
   N'(t) = (dydt[6],dydt[7],dydt[8])
   B'(t) = (dydt[9],dydt[10],dydt[11]). 
   
   The void *params is cast to type "struct intrinsic_polyline_data". 
   
   Note that we borrowed this approach from Gray's 
   Modern Differential Geometry (p. 221).                */
  
{
  double T[3],N[3],B[3],*gp,*Tp,*Np,*Bp;
  double kap,tau;
  int i;
  ktdata_type *ktdata;
  
  /* First, we rename the input and output variables more pleasantly. */

  T[0] = y[3]; T[1] = y[4];  T[2] = y[5];
  N[0] = y[6]; N[1] = y[7];  N[2] = y[8];
  B[0] = y[9]; B[1] = y[10]; B[2] = y[11];

  gp = &(dydt[0]); Tp = &(dydt[3]); Np = &(dydt[6]); Bp = &(dydt[9]);
  
  /* Now we compute kappa and tau. */

  ktdata = (ktdata_type *) params;
  
  kap = ktdata->curvature(t,ktdata->curvedata);
  tau = ktdata->torsion(t,ktdata->curvedata);

  /* And finally, we assign the output results: these are just the */
  /* Frenet equations. */

  for (i=0;i<3;i++) {

    gp[i] = T[i];

    Tp[i] =           kap*N[i];
    Np[i] = -kap*T[i]            + tau*B[i];
    Bp[i] =          -tau*N[i];

  }

  return GSL_SUCCESS;

}

plCurve *ktcurve(plc_vector x0,plc_vector t0,plc_vector n0,
		 int verts,
		 double startt, double endt,
		 ktfunc kappa, ktfunc tau,
		 void *ktdata)

{
  const gsl_odeiv_step_type * T                /* Standard Runge-Kutta 4/5 integration. */
    = gsl_odeiv_step_rkf45;
  
  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc (T, 12);            /* A 12-dimensional stepper of type T */
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new (1e-6, 0.0);     /* Absolute error of 1e-6 */
  gsl_odeiv_evolve * e 
    = gsl_odeiv_evolve_alloc (12);             /* A 12-dimensional step controller */

  /* We also include the spline parameters. */

  int i;
  ktdata_type ktdata_passin;
  gsl_odeiv_system sys = {kt_integration_func, NULL, 12, &ktdata_passin};

  double t = startt;
  double h = 1e-6;
  double y[12];

  plCurve *C;
  plc_vector b0;

  /* First we assign ktdata */

  ktdata_passin.curvature = kappa;
  ktdata_passin.torsion   = tau;
  ktdata_passin.curvedata = ktdata; 

  /* and the initial conditions. */
  
  if (fabs(plc_dot_prod(t0,n0)) > 1e-8*plc_norm(t0)*plc_norm(n0)) {

    kterrcode = 129;
    sprintf(kterrstring,"ktcurve: t0 and n0 are not normal to one another");
    return NULL;

  }

  bool tok,nok;

  t0 = plc_normalize_vect(t0,&tok); n0 = plc_normalize_vect(n0,&nok);

  if (!tok || !nok) { 

    kterrcode = 199;
    sprintf(kterrstring,"t0 or n0 is a zero vector");
    return NULL;

  }
  
  b0 = plc_cross_prod(t0,n0);

  y[0] = x0.c[0]; y[1]  = x0.c[1]; y[2]  = x0.c[2];
  y[3] = t0.c[0]; y[4]  = t0.c[1]; y[5]  = t0.c[2];
  y[6] = n0.c[0]; y[7]  = n0.c[1]; y[8]  = n0.c[2];
  y[9] = b0.c[0]; y[10] = b0.c[1]; y[11] = b0.c[2];

  /* Now we allocate the polyline. */

  int  nv = verts;
  bool open = (1==1);
  int  cc = 1;

  C = plc_new(1,&nv,&open,&cc);

  /* And finally, we loop the stepper. */
  
  C->cp[0].vt[0] = x0;
  
  for (i = 1; i < C->cp[0].nv; i++) {
    
    double ti = (i * (endt) / (C->cp[0].nv-1));
    
    while (t < ti) {
      
      int status = gsl_odeiv_evolve_apply (e, c, s, 
					   &sys, 
					   &t, ti, &h,
					   y);

      if (status != GSL_SUCCESS) {

	kterrcode = 182;
	sprintf(kterrstring,"Ode solver failure %s.\n",gsl_strerror(status));
	
	plc_free(C);
	return NULL;

      }
	
    }
    
    C->cp[0].vt[i].c[0] = y[0]; C->cp[0].vt[i].c[1] = y[1]; C->cp[0].vt[i].c[2] = y[2];

  }
    
  /* We have now created the curve. Free the ODE solver stuff and quit. */

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
 
  return C;

 }
