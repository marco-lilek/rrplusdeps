/*
 *  mapmmodel.c

 This file, part of writhecomp, includes standard C doubles as an arithmetic
 model for writhe computations. We will also probably include a "long double"
 model, to see what effect this has on the computations. 

 */

#include "writhecomp.h"

int MAPM_running = {FALSE};

void init_mapm_model(int precision) 

     /* Procedure just initializes the standard constants PI and np_MINUS. */

{
  M_APM minus;
  M_APM eps;
  M_APM pi;
  M_APM one;
  M_APM four;
  M_APM quarterpi;
  M_APM zero;

  char  eps_string[128];
  char  debug_buf[1024];

  /* We test to see whether the library has already been brought up and shut down. */

  if (MAPM_running) {

    m_apm_trim_mem_usage();

  } else {

    MAPM_running = TRUE;

  }
  
  /* The library MM_PI seems to be failing, so we compute PI for ourselves in the */
  /* requested precision. */

  one  = m_apm_init();
  four = m_apm_init();
  quarterpi = m_apm_init();
  pi = m_apm_init();

  m_apm_set_string(one,"1.0");
  m_apm_set_string(four,"4.0");

  m_apm_arctan(quarterpi,precision+6,one);
  m_apm_multiply(pi,quarterpi,four);

  m_apm_round(pi,precision,pi);

  np_PI = (NUMPTR)(pi);
  m_apm_to_string(debug_buf,25,pi);

  /* The constant np_MINUS will be a little harder. */

  minus = m_apm_init();
  m_apm_set_string(minus,"-1.0");
  
  np_MINUS = (NUMPTR)(minus);
  m_apm_to_string(debug_buf,25,minus);

  /* And np_EPS will be harder still. */

  eps = m_apm_init();
  sprintf(eps_string,"1e-%d",DEC_PRECISION-5);
  m_apm_set_string(eps,eps_string);

  np_EPS = (NUMPTR)(eps);
  m_apm_to_string(debug_buf,25,eps);

  /* We also have np_ZERO. */

  zero = m_apm_init();
  m_apm_set_string(zero,"0.0");
  np_ZERO = (NUMPTR)(zero);

  m_apm_free(one);
  m_apm_free(four);
  m_apm_free(quarterpi);
}

void kill_mapm_model()

     /* Procedure eliminates the standard constants. */

{
 
  kill_mapm(np_MINUS);
  kill_mapm(np_PI);
  kill_mapm(np_EPS);
  kill_mapm(np_ZERO);

  m_apm_free_all_mem();

}

void kill_mapm(NUMPTR a)
     
     /* Procedure frees the memory associated to a. */

{
  m_apm_free((M_APM)(a));
}

NUMPTR double_to_mapm(double x)

     /* Procedure converts a double to a NUMPTR. */

{
  M_APM ret;

  ret = m_apm_init();
  m_apm_set_double(ret,x);

  return (NUMPTR)(ret);
}

double mapm_to_double(NUMPTR x)

     /* Procedure converts a NUMPTR to a mapm. */

{
  char stringval[50];
  double result = {0.0};

  m_apm_to_string(stringval,15,(M_APM)(x));
  
  if(sscanf(stringval,"%lf",&result) != 1) {

    fprintf(stderr,"mapm_to_double: Failed to understand MAPM string %s.\n",stringval);
    exit(1);

  }

  return result;
}

char *mapm_to_string(NUMPTR a)

     /* Procedure converts a mapm number to a string with DEC_PRECISION digits. */

{
  char *output;

  output = calloc(1024,sizeof(char));
  m_apm_to_string(output,DEC_PRECISION,(M_APM)(a));

  return output;
}

NUMPTR mapm_add(NUMPTR a,NUMPTR b)

     /* Procedure adds mapm numbers and rounds the result to DEC_PRECISION */
     /* decimal places. */

{
  M_APM result, scratch;

  result = m_apm_init();
  scratch = m_apm_init();

  m_apm_add(scratch,(M_APM)(a),(M_APM)(b));
  m_apm_round(result,DEC_PRECISION,scratch);
  
  m_apm_free(scratch);
  return result;
}


NUMPTR mapm_mul(NUMPTR a,NUMPTR b)

     /* Procedure multiplies mapm numbers and rounds the result to DEC_PRECISION */
     /* decimal places. */

{
  M_APM result, scratch;

  result = m_apm_init();
  scratch = m_apm_init();

  m_apm_multiply(scratch,(M_APM)(a),(M_APM)(b));
  m_apm_round(result,DEC_PRECISION,scratch);
  
  m_apm_free(scratch);
  return result;
}


NUMPTR mapm_div(NUMPTR a,NUMPTR b)

     /* Procedure divides mapm numbers and rounds the result to DEC_PRECISION */
     /* decimal places. */

{
  M_APM result, scratch;

  result = m_apm_init();
  scratch = m_apm_init();

  m_apm_divide(scratch,DEC_PRECISION,(M_APM)(a),(M_APM)(b));
  m_apm_round(result,DEC_PRECISION,scratch);
  
  m_apm_free(scratch);
  return result;
}

NUMPTR mapm_acos(NUMPTR theta)

     /* Procedure takes the arc cosine of theta using the mapm
	library function. In principle, this is accurate to the 
	number of digits in DEC_PRECISION. */
{
  M_APM result;

  result = m_apm_init();
  m_apm_arccos(result,DEC_PRECISION,(M_APM)(theta));

  return (NUMPTR)(result);
}


NUMPTR mapm_sqrt(NUMPTR x)

     /* Procedure takes the square root of x using the mapm library function. 
	This should be accurate to within DEC_PRECISION digits. */

{
  M_APM result;

  result = m_apm_init();
  m_apm_sqrt(result,DEC_PRECISION,(M_APM)(x));

  return (NUMPTR)(result);
}

int mapm_leq(NUMPTR a, NUMPTR b) 

     /* Procedure compares doubles. */

{
  if (m_apm_compare((M_APM)(a),(M_APM)(b)) < 1) {

    return TRUE;

  } else {

    return FALSE;

  }
}

