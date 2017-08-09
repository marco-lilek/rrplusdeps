/* 

   kthelix.c : This test program uses the ktcurve library to construct a helix.

*/

#include "plCurve.h"
#include "ktcurve.h"

#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#ifdef HAVE_STDLIB_H
  #include "stdlib.h"
#endif

double kappatest(double t, void *data)
{
  return 2;
}

double tautest(double t, void *data)
{
  return 1;
}

int main()
{
  plCurve *newhelix;
  plc_vector x0 = {{1,0,0}}, t0 = {{0,1,0}}, n0 = {{-1,0,0}};

  printf("kthelix: Generating a 300 vertex helix from curvature and torsion.\n");
  
  newhelix = ktcurve(x0,t0,n0,300,0,15,kappatest,tautest,NULL);

  if (newhelix == NULL) {

    printf("kthelix: Error: %s.\n",kterrstring);
    exit(1);

  }

  FILE *outfile;

  outfile = fopen("helix.vect","w");
  plc_write(outfile,newhelix);

  fclose(outfile);
  plc_free(newhelix);

  printf("kthelix: Curve generated.\n");
  exit(0);
}
  
