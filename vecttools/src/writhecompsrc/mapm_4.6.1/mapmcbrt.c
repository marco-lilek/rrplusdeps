
/* 
 *  M_APM  -  mapmcbrt.c
 *
 *  Copyright (C) 2000 - 2002   Michael C. Ring
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
 *      $Id: mapmcbrt.c,v 1.5 2002/11/03 21:34:34 mike Exp $
 *
 *      This file contains the CBRT (cube root) function.
 *
 *      $Log: mapmcbrt.c,v $
 *      Revision 1.5  2002/11/03 21:34:34  mike
 *      Updated function parameters to use the modern style
 *
 *      Revision 1.4  2000/10/30 16:42:22  mike
 *      minor speed optimization
 *
 *      Revision 1.3  2000/07/11 18:03:39  mike
 *      make better estimate for initial precision
 *
 *      Revision 1.2  2000/04/08 18:34:35  mike
 *      added some more comments
 *
 *      Revision 1.1  2000/04/03 17:58:04  mike
 *      Initial revision
 */

#include "m_apm_lc.h"

/****************************************************************************/
void	m_apm_cbrt(M_APM rr, int places, M_APM aa)
{
M_APM   last_x, guess, tmpN, tmp7, tmp8, tmp9;
int	ii, nexp, bflag, tolerance, maxp, local_precision;

/* result is 0 if input is 0 */

if (aa->m_apm_sign == 0)
  {
   M_set_to_zero(rr);
   return;
  }

last_x = M_get_stack_var();
guess  = M_get_stack_var();
tmpN   = M_get_stack_var();
tmp7   = M_get_stack_var();
tmp8   = M_get_stack_var();
tmp9   = M_get_stack_var();

/* compute the cube root of the positive number, we'll fix the sign later */

m_apm_absolute_value(tmpN, aa);
nexp = 0;

/* 
    normalize the input number (make the exponent near 0) so
    the 'guess' function will not over/under flow on large
    magnitude exponents.
*/

if (aa->m_apm_exponent > 2)
  {
   nexp = aa->m_apm_exponent / 3;
   tmpN->m_apm_exponent -= 3 * nexp;
  }

if (aa->m_apm_exponent < -2)
  {
   nexp = aa->m_apm_exponent / 3 - 1;
   tmpN->m_apm_exponent -= 3 * nexp;
  }

M_get_cbrt_guess(guess, tmpN);

tolerance = places + 4;
maxp      = places + 16;
bflag     = FALSE;

m_apm_multiply(tmp9, guess, guess);
m_apm_multiply(tmp8, guess, tmp9);
m_apm_subtract(tmp9, tmpN, tmp8);

if ((local_precision = 16 - tmp9->m_apm_exponent) < 30)
  local_precision = 30;

m_apm_negate(last_x, MM_Ten);

/*   Use the following iteration to calculate the cube root :

                   1               N
         X     =  --- * [ 2 * X + --- ]
          n+1      3              X^2
*/

ii = 0;

while (TRUE)
  {
   m_apm_multiply(tmp8, MM_Two, guess);
   m_apm_multiply(tmp7, guess, guess);
   m_apm_divide(tmp9, local_precision, tmpN, tmp7);
   m_apm_add(tmp7, tmp8, tmp9);
   m_apm_divide(guess, local_precision, tmp7, MM_Three);

   if (bflag)
     break;

   /* force at least 2 iterations so 'last_x' has valid data */

   if (ii != 0)
     {
      m_apm_subtract(tmp8, guess, last_x);

      if (tmp8->m_apm_sign == 0)
        break;

      if ((-4 * tmp8->m_apm_exponent) > tolerance)
        bflag = TRUE;
     }

   if (ii == 0)
     local_precision *= 2;
   else
     local_precision += 2 - 2 * tmp8->m_apm_exponent;

   if (local_precision > maxp)
     local_precision = maxp;
  
   m_apm_copy(last_x, guess);
   ii++;
  }

m_apm_round(rr, places, guess);
rr->m_apm_exponent += nexp;
rr->m_apm_sign = aa->m_apm_sign;
M_restore_stack(6);
}
/****************************************************************************/
