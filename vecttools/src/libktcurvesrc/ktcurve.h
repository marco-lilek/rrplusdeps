/* ktcurve.h.  Generated by configure.  */
/*
 *
 * A library for generating plCurves from curvature and torsion using GSL.
 *
 *  $Id: ktcurve.h,v 1.1 2008/01/24 05:05:01 cantarel Exp $
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

#ifndef KTCURVE_H
#define KTCURVE_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif
  
  #include "plCurve.h"

  extern int  kterrcode;
  extern char kterrstring[1024];

  typedef double ktfunc(double, void *);

  /* This is the general prototype for the curvature and torsion
     functions.  These functions will each be passed a parameter value
     between startt and endt and the void *ktdata.

     They must evaluate to a value for curvature and torsion. The
     solver assumes that these functions are at least continuous in
     the parameter, but it will do the best it can with discontinuous
     k and t functions. */
 
  plCurve *ktcurve(plc_vector x0,plc_vector t0,plc_vector n0,
		   int verts,
		   double startt, double endt,
		   ktfunc kappa, ktfunc tau,
		   void *ktdata);
  
  /* Integrates the curvature and torsion functions kappa and tau over
     (startt,endt) to create a plCurve with verts vertices. The curve
     starts at x0 with initial tangent vector t0 and initial normal
     vector n0. The void ptr ktdata is passed to kappa and tau
     whenever they are evaluated. 

     If the integration fails (usually due to problems with kappa or
     tau functions), the procedure returns NULL and sets kterrcode and
     kterrstring. */

#if (__cplusplus || c_plusplus)
};
#endif
#endif /* KTCURVE_H */