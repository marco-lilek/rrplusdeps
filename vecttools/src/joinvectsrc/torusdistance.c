#include"joinvect.h"
#include<nplCurve.h>
#include<math.h>

/* 

  torusdistance gives the distance between points in the flat torus of the dimension
  (call it N) of the given nplc_vectors. This is updated to be quite fast.

*/

/* Finds the longest edge in an nplCurve */
void longest_edge(const nplCurve * const L, int *cp, int *vt0, int *vt1, double *len,
		  double twidth)

{
  int i,j;
  double thisd;

  *len = -1.0;

  for(i=0;i<L->nc;i++) {

    for(j=0;j<L->cp[i].nv;j++) {

      thisd = torus_distance(L->cp[i].vt[j],L->cp[i].vt[j+1],twidth);

      if (thisd > *len) {

	*cp = i; *vt0 = j; *vt1 = (j+1) % L->cp[i].nv; *len = thisd;

      }

    }

  }

}

double cdist(double twidth,double a, double b)
{
  
  double cd;

  cd = fabs(a - b);
  if (cd > fabs(a - (b + twidth))) { cd = fabs(a - (b + twidth)); }
  if (cd > fabs(a - (b - twidth))) { cd = fabs(a - (b - twidth)); }

  return cd;

}

double torus_distance(nplc_vector A, nplc_vector B,double twidth)
{
  double pd=0,cd;
  int i;

  for(i=0;i<A.n;i++) { 

    cd = cdist(twidth,A.c[i],B.c[i]);
    pd += cd * cd;
  
  }

  pd = sqrt(pd);
  
  return pd;

}

