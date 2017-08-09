#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#include "tube.h"
#include "plCurve.h"

plc_vector tube_tangent(plCurve *L,int cmp,int vert,bool *ok)

/* Finds the tangent to link->cp[cp].vt as the average of the tangents on either side, NOT COUNTING LENGTHS */

{
  if (L->cp[cmp].open) {
    /* For open strands, take the only tangent we have. */
    if (vert == 0) {
      return plc_normalize_vect(
          plc_vect_diff(L->cp[cmp].vt[vert+1],L->cp[cmp].vt[vert]),
          ok);
    } else if (vert == L->cp[cmp].nv-1) {
      return plc_normalize_vect(
          plc_vect_diff(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert-1]),
          ok);
    }
  }

  /* We now know that either we are on a closed
     component, or we are not at an endpoint.   */

  plc_vector ftan, rtan;
  bool isok;

  ftan = plc_normalize_vect(plc_vect_diff(L->cp[cmp].vt[vert+1],L->cp[cmp].vt[vert]),&isok);
  rtan = plc_normalize_vect(plc_vect_diff(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert-1]),&isok);

  return plc_normalize_vect(plc_vect_sum(ftan,rtan),ok);

} /* tube_tangent */
  

plCurve *split_sharp_corners(plCurve *L)
  
/* At each turning angle above 5 degrees, split the incident edges in half. */
  
{
  int *nv,*cc;
  bool *open,*vertcolors;
  double maxturn = 3.1415926 * 5.0 / 180.0;
  int cp,vt;
  plCurve *newL;

  /* Create a copy with room for new vertices in place */

  nv = calloc(L->nc,sizeof(int));
  cc = calloc(L->nc,sizeof(int));
  open = calloc(L->nc,sizeof(bool));
  vertcolors = calloc(L->nc,sizeof(bool));

  for(cp=0;cp<L->nc;cp++) {

    nv[cp] = L->cp[cp].nv * 2;
    vertcolors[cp] = (L->cp[cp].cc == L->cp[cp].nv);
    cc[cp] = vertcolors[cp] ? 2*L->cp[cp].nv : L->cp[cp].cc;
    open[cp] = L->cp[cp].open;
    

  }

  newL = plc_new(L->nc,nv,open,cc);

  free(nv);
  free(cc);
  free(open);

  /* Now go through the list, adding vertices as needed */

  for(cp=0;cp<L->nc;cp++) {

    /* We are looping over edges here, and adding the trailing vertex from each edge to L on each step. */

    newL->cp[cp].nv = 0; // We will add vertices to this buffer one-by-one.
    if (vertcolors[cp]) { newL->cp[cp].cc = 0; }
    else if (newL->cp[cp].cc > 0) { newL->cp[cp].clr[0] = L->cp[cp].clr[0]; }

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      newL->cp[cp].vt[newL->cp[cp].nv++] = L->cp[cp].vt[vt]; // Add the trailing vertex
      if (vertcolors[cp]) { newL->cp[cp].clr[newL->cp[cp].cc++] =  L->cp[cp].clr[vt]; } // and it's color (if assigned)

      if (plc_turning_angle(L,cp,vt,NULL) > maxturn || plc_turning_angle(L,cp,vt+1,NULL) > maxturn) {

	newL->cp[cp].vt[newL->cp[cp].nv++] = plc_vweighted(0.5,L->cp[cp].vt[vt],L->cp[cp].vt[vt+1]); // Add midpoint
	if (vertcolors[cp]) { newL->cp[cp].clr[newL->cp[cp].cc++] =  L->cp[cp].clr[vt]; } // and old color (if assigned)

      }
	
    }

  }
	
  plc_fix_wrap(newL);
  free(vertcolors);
    
  return newL;
    
}




plCurve *plCurve_equalize_sides(plCurve *L, int *target_verts)

  /* Procedure divides each polyline into target_verts sides, all
     of nearly equal length, by replacing the vertices of the ith pline
     with target_verts[i] vertices, spaced equally _along the curve_. 

     The procedure generates a new link with the target_vert counts.  
     We lose the color information from the original link in this version. */

{
  int             i,j,cmp;
  plc_vector  disp;
  double          side_length,length_so_far,t,seg_length;
  double         *component_lengths;

  plCurve   *new_L;
  int       *cc;
  bool      *open;

  assert(L != NULL);
  assert(target_verts != NULL);

  /* Allocate new memory for the revised link. */

  open = (bool *)(calloc(L->nc,sizeof(bool)));
  cc   = (int *)(calloc(L->nc,sizeof(int)));

  for(i=0;i<L->nc;i++) {

   open[i] = L->cp[i].open;
   cc[i]   = L->cp[i].cc;

  }

  new_L = plc_new(L->nc,target_verts,open,cc);
  free(open);
  free(cc);

  /* Now start computing the new link */

  component_lengths = (double *)(calloc(L->nc,sizeof(double)));
  plc_arclength(L,component_lengths);

  for(cmp=0;cmp<L->nc;cmp++) {

    side_length = component_lengths[cmp]/
      (double)(L->cp[cmp].open ? target_verts[cmp] - 1: target_verts[cmp]);

    new_L->cp[cmp].vt[0] = L->cp[cmp].vt[0];
    length_so_far = 0; t = 0; j = 1;

    /* We are now ready to loop through pline. */

    for(i=0;j<target_verts[cmp];i++) {
      
      disp = plc_vect_diff(L->cp[cmp].vt[i+1],L->cp[cmp].vt[i]);
      seg_length = plc_M_norm(disp);
      t = 0;
      
      while (t < 1) { /* While on the current segment of pline */
    
	if (side_length - length_so_far < (1-t)*seg_length) {
	  
	  /* If there is enough length here to end the new segment... */
	  
	  t += (side_length - length_so_far)/(seg_length);
	  plc_M_vlincomb(new_L->cp[cmp].vt[j],1,L->cp[cmp].vt[i],t,disp);
	  j++; length_so_far = 0;
	  
	} else {
	  
	  /* If not, we'll advance length_so_far, reset t, and increment i */

	  length_so_far += (1-t)*seg_length;
	  t = 1;
	  
	}
	
      }

    }
    
    /* Now, we should have generated all the new vertices. */
    /* If the polyline is.open, we need to set the last vert manually. */

    if (L->cp[cmp].open) {
    
      new_L->cp[cmp].vt[target_verts[cmp]-1] = L->cp[cmp].vt[L->cp[cmp].nv-1];

    }
    
  }

  plc_fix_wrap(new_L);
  return new_L;

}


plc_vector *random_framezeros(plCurve *link)

/* Procedure generates random initial vectors for orthonormal
   frames on the components of link. */

{
  plc_vector *fz_buf;
  plc_vector T;
  int i,j;
  double dprod;
  
  fz_buf = calloc(link->nc,sizeof(plc_vector));
  assert(fz_buf != NULL);

  srand((int)(time(NULL)));

  bool ok;
  
  for(i=0;i<link->nc;i++) {

    T = tube_tangent(link,i,0,&ok);
    
    if (!ok) {

      printf("tube: Problem generating tangent for vert 0 of component %d of link.\n"
	     "      This could be a plCurve problem. Suggest you recompile plCurve\n"
             "      and tube.\n",i);
      exit(1);

    }

    for(j=0;j<200;j++) {/* We worry about picking R and fz_buf[i] colinear. */

      fz_buf[i] = plc_random_vect();
      dprod = plc_M_dot(fz_buf[i],T);
      if (fabs(dprod) < 0.9) break;

    }
      
    assert(j < 200); /* And assume we've succeeded after 200 tries. */
    
    plc_M_scale_vect(dprod,T);		/* Ensure that dot(T,fz_buf[i]) = 0. */
    plc_M_sub_vect(fz_buf[i],T);
    fz_buf[i] = plc_normalize_vect(fz_buf[i],&ok);

    assert(ok);
  }

  return fz_buf;

}

void plCurve_bishop_frame(plCurve *link,
			       plc_vector *framezero,
			       plCurve **frameA, plCurve **frameB)

/* Procedure creates the "Bishop frame" on each component of the link
   with initial vectors given by "framezeros". */

{
  plc_vector T;
  int cp, vt;
  double dprod;

  /* First, we allocate space for the frame. */

  (*frameA) = plc_copy(link);
  (*frameB) = plc_copy(link);

  /* Now we are ready to fill it in. */

  bool ok;
			       
  for(cp=0;cp < link->nc;cp++) {

    T = tube_tangent(link,cp,0,&ok);
    assert(ok);

    if (plc_M_dot(T,framezero[cp]) > 1e-6) {

      /* The initial vector is not normal to T,
         so we fix it. */

      dprod = plc_M_dot(T,framezero[cp]);
      plc_M_vmadd(framezero[cp],-dprod,T);
      framezero[cp] = plc_normalize_vect(framezero[cp],&ok);
      assert(ok);

    }

    (*frameA)->cp[cp].vt[0] = framezero[cp];
    (*frameB)->cp[cp].vt[0] = plc_cross_prod(T,framezero[cp]);
    (*frameB)->cp[cp].vt[0] = plc_normalize_vect((*frameB)->cp[cp].vt[0],&ok);
    assert(ok);
    
    for(vt=1;vt < link->cp[cp].nv;vt++) {

      T = tube_tangent(link,cp,vt,&ok);
      assert(ok);

      (*frameA)->cp[cp].vt[vt] = (*frameA)->cp[cp].vt[vt-1];

      dprod = plc_M_dot(T,(*frameA)->cp[cp].vt[vt]);
      plc_M_vmadd((*frameA)->cp[cp].vt[vt],-dprod,T);
      (*frameA)->cp[cp].vt[vt] = plc_normalize_vect((*frameA)->cp[cp].vt[vt],
						    &ok);
      assert(ok);

      (*frameB)->cp[cp].vt[vt] = 
	plc_normalize_vect(plc_cross_prod(T,(*frameA)->cp[cp].vt[vt]),&ok);
      assert(ok);

    }
      
  }

  /* We've now generated the frame. */

}

void plc_force_closed( plCurve *inLink )

     /* Forces each open component of inLink to be closed. */

{
  int cItr, fixItr;
  plc_vector error;
  double efrac;

  for(cItr=0;cItr<inLink->nc;cItr++) {

    if (inLink->cp[cItr].open) { 

      error = plc_vect_diff(inLink->cp[cItr].vt[inLink->cp[cItr].nv-1],
			    inLink->cp[cItr].vt[0]);

      for( fixItr=0; fixItr < inLink->cp[cItr].nv; fixItr++) {

	efrac = (double)fixItr/(double)inLink->cp[cItr].nv-1;

	inLink->cp[cItr].vt[fixItr] = plc_vmadd(inLink->cp[cItr].vt[fixItr],
						efrac,error);

      }

      /* At this point, the last and first vertices are the same. */

      inLink->cp[cItr].nv--;
      inLink->cp[cItr].open = false;
      
    }

  }

  plc_fix_wrap(inLink);

}
						
      
      
      

      

      
   
