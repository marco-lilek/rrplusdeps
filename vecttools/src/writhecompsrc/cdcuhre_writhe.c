#include <cubature.h>
#include <plCurve.h>
#include <math.h>
#include <stdlib.h>

/* 

   Numerical integration of writhe using the cubature recursive subdivision algorithm.

*/

struct component_pair {
  
  int cp[2];
  plCurve *L;

};

struct spline_component_pair {
  
  int cp[2];
  plc_spline *spL;

};

/* Denotes a pair of components to integrate over. */

void writhe_integrand(unsigned dim,const double *pt,void *data,unsigned nfuncs,double *val)

/* Evaluates the writhe integrand at coordinates (pt[0],pt[1]) on the curve given by data. */
/* The integration is NOT with respect to arclength, so we are integrating

   (T(t0) x T(t1) . (L(t1) - L(t0)))/|L(t1) - L(t0)|^3 * |L'(t0)| |L'(t1)| dt0 dt1

   We can use trilinearity of the triple product to replace the T(t0) |L'(t0)| and T(t1) |L'(t1)| 
   terms with L'(t0) and L'(t1). Now we are going to let t range from 0 to the number of vertices
   in the curve, so that each edge is parameter length 1. 

   This means that L'(t0) is always just the edge vector of the curve. 
   We assume that we are given a pair of components and a plCurve in data.

*/

{
  plc_vector tan[2],pos[2];
  plc_vector diff;
  int edges[2];
  struct component_pair *cPair;
  int i;

  cPair = (struct component_pair *)(data);

  for(i=0;i<2;i++) {
    
    edges[i] = floor(pt[i]);
    tan[i] = plc_vect_diff(cPair->L->cp[cPair->cp[i]].vt[edges[i]+1],
			   cPair->L->cp[cPair->cp[i]].vt[edges[i]]);
    
    pos[i] = plc_vweighted(pt[i] - edges[i],
			   cPair->L->cp[cPair->cp[i]].vt[edges[i]],
			   cPair->L->cp[cPair->cp[i]].vt[edges[i]+1]);
  }

  if (edges[0] == edges[1] && cPair->cp[0] == cPair->cp[1]) { /* on the same edge */ 

    val[0] = 0;

  } else {
			     
    diff = plc_vect_diff(pos[0],pos[1]);
    diff = plc_scale_vect(1.0/pow(plc_norm(diff),3.0),diff);

    val[0] = plc_dot_prod(plc_cross_prod(tan[0],tan[1]),diff); 
  
  }

}

double fast_writhe(plCurve *L)

/* Computes writhe of L with adaptive subdivision. */

{
  int i,j;
  double writhe_accum = 0.0,pair_writhe,pair_err;
  struct component_pair cpData;
  double a[2] = {0,0},b[2];
  double oneover4pi = 0.07957747154594766788444188;
  int *edges;

  /* We first do a little setup */

  edges = calloc(L->nc,sizeof(int));
  plc_edges(L,edges);
  cpData.L = L;

  /* Now we are ready to compute. */

  for(i=0;i<L->nc;i++) {

    for(j=0;j<=i;j++) {

      cpData.cp[0] = i; cpData.cp[1] = j;
      b[0] = edges[i]; b[1] = edges[j];
      
      adapt_integrate(1,writhe_integrand,(void *)&cpData,2,a,b,0,1e-4,1e-4,&pair_writhe,&pair_err);
      pair_writhe *= oneover4pi;

      printf("Computed pairwise writhe (link) for cpts %d <-> %d of %g with error %g.\n",
	     i,j,pair_writhe,pair_err);

      writhe_accum += ((i == j) ? 1 : 2) * pair_writhe;

    }

  }

  printf("Total writhe: %g\n\n",writhe_accum);
  free(edges);

  return writhe_accum;
}

      
  
void writhe_spline_integrand(unsigned dim,const double *pt,void *data,unsigned nfuncs,double *val)

/* Evaluates the writhe integrand at coordinates (pt[0],pt[1]) on the spline curve given by data. */
/* The integration is NOT with respect to arclength, so we are integrating

   (T(t0) x T(t1) . (L(t1) - L(t0)))/|L(t1) - L(t0)|^3 * |L'(t0)| |L'(t1)| dt0 dt1

   We can use trilinearity of the triple product to replace the T(t0) |L'(t0)| and T(t1) |L'(t1)| 
   terms with L'(t0) and L'(t1). Now we are going to let t range from 0 to the number of vertices
   in the curve, so that each edge is parameter length 1. 

   This means that L'(t0) is always just the edge vector of the curve. 
   We assume that we are given a pair of components and a plCurve in data.

*/

{
  plc_vector tan[2],pos[2];
  plc_vector diff;
  struct spline_component_pair *cPair;
  int i;

  cPair = (struct spline_component_pair *)(data);

  for(i=0;i<2;i++) {

    pos[i] = plc_sample_spline(cPair->spL,cPair->cp[i],pt[i]);
    tan[i] = plc_spline_tangent(cPair->spL,cPair->cp[i],pt[i]);

  }

  if (fabs(pt[0] - pt[1]) < 1e-8 && cPair->cp[0] == cPair->cp[1]) { /* too close */ 

    val[0] = 0;

  } else {
			     
    diff = plc_vect_diff(pos[0],pos[1]);
    diff = plc_scale_vect(1.0/pow(plc_norm(diff),3.0),diff);

    val[0] = plc_dot_prod(plc_cross_prod(tan[0],tan[1]),diff); 
  
  }

}

double fast_spline_writhe(plCurve *L)

/* Computes writhe of L with adaptive subdivision. */

{
  int i,j;
  double writhe_accum = 0.0,pair_writhe,pair_err;
  struct spline_component_pair cpData;
  double a[2] = {0,0},b[2];
  double oneover4pi = 0.07957747154594766788444188;
  double *arclength;
  bool ok;

  /* We first do a little setup */

  arclength = calloc(L->nc,sizeof(double));
  plc_arclength(L,arclength);
  cpData.spL = plc_convert_to_spline(L,&ok);

  if (!ok) {

    fprintf(stderr,"writhecomp: Died attempting to convert to spline.\n");
    exit(1);

  }

  /* Now we are ready to compute. */

  for(i=0;i<L->nc;i++) {

    for(j=0;j<=i;j++) {

      cpData.cp[0] = i; cpData.cp[1] = j;
      b[0] = arclength[i]; b[1] = arclength[j];
      
      adapt_integrate(1,writhe_spline_integrand,(void *)&cpData,2,a,b,0,1e-4,1e-4,&pair_writhe,&pair_err);
      pair_writhe *= oneover4pi;

      printf("Computed pairwise writhe (link) for cpts %d <-> %d of %g with error %g by spline method.\n",
	     i,j,pair_writhe,pair_err);

      writhe_accum += ((i == j) ? 1 : 2) * pair_writhe;

    }

  }

  plc_spline_free(cpData.spL);
  printf("Total writhe: %g\n\n",writhe_accum);
  free(arclength);

  return writhe_accum;
}
