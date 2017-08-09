#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include"maxmin.h"

int intmax(int nargs, ... )

/* intmax finds the max of an arbitrary collection of <nargs> integer 
   arguments. */

{
  va_list ap;
  int     themax,thisarg,i;

  assert(nargs >= 1);
  va_start(ap, nargs);

  themax = va_arg(ap,int);

  for(i=1;i<nargs;i++) {

    thisarg = va_arg(ap,int);
    themax = (themax > thisarg ? themax : thisarg);

  }

  va_end(ap);

  return themax;
}

int intmin(int nargs, ... )

/* intmin finds the min of an arbitrary collection of <nargs> integer 
   arguments. */

{
  va_list ap;
  int     themin,thisarg,i;

  assert(nargs >= 1);
  va_start(ap, nargs);

  themin = va_arg(ap,int);

  for(i=1;i<nargs;i++) {

    thisarg = va_arg(ap,int);
    themin = (themin < thisarg ? themin : thisarg);

  }

  va_end(ap);

  return themin;
}

double doublemax(int nargs, ... )

/* doublemax finds the max of an arbitrary collection of <nargs> double
   arguments. */

{
  va_list ap;
  double     themax,thisarg,i;

  assert(nargs >= 1);
  va_start(ap, nargs);

  themax = va_arg(ap,double);

  for(i=1;i<nargs;i++) {

    thisarg = va_arg(ap,double);
    themax = (themax > thisarg ? themax : thisarg);

  }

  va_end(ap);

  return themax;
}

double doublemin(int nargs, ... )

/* doublemin finds the min of an arbitrary collection of <nargs> integer 
   arguments. */

{
  va_list ap;
  double     themin,thisarg,i;

  assert(nargs >= 1);
  va_start(ap, nargs);

  themin = va_arg(ap,double);

  for(i=1;i<nargs;i++) {

    thisarg = va_arg(ap,double);
    themin = (themin < thisarg ? themin : thisarg);

  }

  va_end(ap);

  return themin;
}



  
