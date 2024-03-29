
/* 
 *  M_APM  -  mapm_lg2.c
 *
 *  Copyright (C) 2003   Michael C. Ring
 *
 *  Permission to use, copy, and distribute this software and its
 *  documentation for any purpose with or without fee is hereby granted, 
 *  provided that the above copyright notice appear in all copies and 
 *  that both that copyright notice and this permission notice appear 
 *  in supporting documentation.
 *
 *  Permission to modify the software is granted, but not the right to
 *  distribute the modified code.  Modifications are to be distributed 
 *  as patches to released version.
 *  
 *  This software is provided "as is" without express or implied warranty.
 */

/*
 *      $Id: mapm_lg2.c,v 1.4 2003/04/09 20:21:29 mike Exp $
 *
 *      This file contains the iterative function to compute the LOG
 *	This is an internal function to the library and is not intended
 *	to be called directly by the user.
 *
 *      $Log: mapm_lg2.c,v $
 *      Revision 1.4  2003/04/09 20:21:29  mike
 *      fix rare corner condition by intentionally inducing a
 *      10 ^ -5 error in the initial guess.
 *
 *      Revision 1.3  2003/03/31 22:13:15  mike
 *      call generic error handling function
 *
 *      Revision 1.2  2003/03/30 21:27:22  mike
 *      add comments
 *
 *      Revision 1.1  2003/03/30 21:18:07  mike
 *      Initial revision
 */

#include "m_apm_lc.h"
#include <math.h>

/****************************************************************************/

/*
 *      compute rr = log(nn)
 *
 *	input is assumed to not exceed the exponent range of a normal
 *	'C' double ( |exponent| must be < 308)
 */

/****************************************************************************/
void	M_log_basic_iteration(M_APM rr, int places, M_APM nn)
{
M_APM   tmp0, tmp1, tmp2, tmp3, guess;
int	ii, maxiter, maxp, tolerance, local_precision;

guess = M_get_stack_var();
tmp0  = M_get_stack_var();
tmp1  = M_get_stack_var();
tmp2  = M_get_stack_var();
tmp3  = M_get_stack_var();

M_get_log_guess(guess, nn);

tolerance       = -(places + 4);
maxp            = places + 16;
local_precision = 18;

/*
 *      compute the maximum number of iterations
 *	that should be needed to calculate to
 *	the desired accuracy.  [ constant below ~= 1 / log(3) ]
 */

maxiter = (int)(log((double)(places + 2)) * 0.91024) + 3;

if (maxiter < 5)
  maxiter = 5;

/*    Use the following iteration to solve for log :

                        exp(X) - N 
      X     =  X - 2 * ------------
       n+1              exp(X) + N 

   
      this is a cubically convergent algorithm 
      (each iteration yields 3X more digits)
*/

ii = 0;

while (TRUE)
  {
   m_apm_exp(tmp1, local_precision, guess);

   m_apm_subtract(tmp3, tmp1, nn);
   m_apm_add(tmp2, tmp1, nn);

   m_apm_divide(tmp1, local_precision, tmp3, tmp2);
   m_apm_multiply(tmp0, MM_Two, tmp1);
   m_apm_subtract(tmp3, guess, tmp0);

   if (ii != 0)
     {
      if (((3 * tmp0->m_apm_exponent) < tolerance) || (tmp0->m_apm_sign == 0))
        break;
     }

   if (++ii == maxiter)
     {
      M_apm_log_error_msg(M_APM_RETURN, 
        "Warning! ... \'M_log_basic_iteration\', max iteration count reached");
      break;
     }

   m_apm_round(guess, local_precision, tmp3);

   local_precision *= 3;

   if (local_precision > maxp)
     local_precision = maxp;
  }

m_apm_round(rr, places, tmp3);
M_restore_stack(5);
}
/****************************************************************************/
