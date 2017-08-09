/*

  tsporder.c : Travelling salesman ordering for curves passed to joinvect.
               Obviously, this is only going to work if ncurves is small. 

*/

#include "joinvect.h"

#define dist(i,j) glob_distances[(i)*ncurves+(j)]

double *glob_distances = 0;

void tsp_recurser(int level,int ncurves,bool *used,
		  double thislen,double *minlen,
		  int order[],int workingorder[])

/* Searches for next links to add at this level. */

{

  int i;

  if (level == ncurves) { /* We've completed a circuit. Check if this wins. */

    if (thislen < *minlen) { 

      for(i=0;i<ncurves;i++) { order[i] = workingorder[i]; }
      *minlen = thislen;

    }

    return;

  }

  /* We are in the middle of building. Loop over the possibilities for this var. */

  for(i=0;i<ncurves;i++) {

    if (!used[i]) { /* if already used, ignore this option */

      if (thislen + dist(level,i) < *minlen) { /* if too long, ignore this option */

	used[i] = true;
	workingorder[level] = i;
	
	tsp_recurser(level+1,ncurves,used,thislen+dist(level,i),minlen,order,workingorder);
	
	used[i] = false;

      }

    }

  }

}

/* void nplc_reverse(nplCurve *L) */

/* /\* Reverses the orientation of L. *\/ */

/* { */
/*   int cp,vt; */
/*   nplc_vector swap; */

/*   for(cp=0;cp<L->nc;cp++) { */

/*     for(vt=0;vt<L->cp[cp].nv/2;vt++) { */

/*       swap = L->cp[cp].vt[vt]; */
/*       L->cp[cp].vt[vt] = L->cp[cp].vt[L->cp[cp].nv-vt-1]; */
/*       L->cp[cp].vt[L->cp[cp].nv-vt-1] = swap; */

/*     } */
/*   } */
/* } */

nplc_vector nplc_start(nplCurve *L,int cp)

{
  return L->cp[cp].vt[0];
}

nplc_vector nplc_end(nplCurve *L,int cp)
{
  return L->cp[cp].vt[L->cp[cp].nv-1];
}

nplc_vector nplc_sf(nplCurve *L,int cp,int sf)
{
  if (sf == 0) {

    return nplc_start(L,cp);

  } else {

    return nplc_end(L,cp);

  }
}

void tsp_order(int ncurves,nplCurve **curves,int *order,double twidth)

/* Orders the curves so that every end vector is joined to a start vector. */
/* Note that we may have to reverse certain curves in order to get the best */
/* result. The resulting array is to be read as follows:

   2 4 1 0 3

   means 0 -> 2, 1 -> 4, 2 -> 1, 3 -> 0, and 4 -> 3, creating loop

   0 -> 2 -> 1 -> 4 -> 3 -> 0.

   Curves that are reversed are reversed in the nplCurve array. */

{
  int i,j,k;
  bool   *used;

  glob_distances = calloc(ncurves*ncurves,sizeof(double));
  used = calloc(ncurves,sizeof(double));

  double minlen = 1e80;
  int    *reversed,*bestreversed,*workingorder;

  reversed = calloc(ncurves,sizeof(int));
  workingorder = calloc(ncurves,sizeof(int));
  bestreversed = calloc(ncurves,sizeof(int));

  for(k=0;k<pow(2,ncurves);k++) {

    /* We rebuild the distance array depending on the orientation of the curves. */

    for(i=0;i<ncurves;i++) {
      
      for(j=0;j<ncurves;j++) {
	
	dist(i,j) = torus_distance(
				 reversed[i] ? nplc_start(curves[i],0) : nplc_end(curves[i],0),
				 reversed[j] ? nplc_end(curves[j],0) : nplc_start(curves[j],0),
				 twidth);
	
      }
      
    }

    for(i=0;i<ncurves;i++) { workingorder[i] = 0; }

    double old_minlen;
    old_minlen = minlen;

    tsp_recurser(0,ncurves,used,0.0,&minlen,order,workingorder);

    if (minlen < old_minlen) { /* We won! Quick, store the reversals... */

      for(i=0;i<ncurves;i++) { bestreversed[i] = reversed[i]; }

    }

    /* Now increment the reversed keys */

    for(j=0;j<ncurves;j++) {

      reversed[j]++;                      // increment this bit
      if (reversed[j] == 1) { break; }    // if no carry, we're done     
      reversed[j] = 0;                    // if carry, we set this bit to 0 and continue to the next one.
      
    }


  }

  /* Now we're done. Actually reverse the needed curves. */

  for(i=0;i<ncurves;i++) {

    if (bestreversed[i]) { nplc_reverse(curves[i]); }

  }

  free(workingorder);
  free(used);
  free(glob_distances);
  glob_distances = NULL;

}

struct endpt {

  int curve;
  int sf;

};

struct distrec {
  
  struct endpt pts[2];
  double dist;
  
};
  
int distrec_cmp(const void *a,const void *b)
{
  const struct distrec *da = (const struct distrec *) a;
  const struct distrec *db = (const struct distrec *) b;
     
  return (da->dist > db->dist) - (da->dist < db->dist);
}

void greedy_order(int ncurves,nplCurve **curves,int *order,double twidth)

{	
  bool           *used[2];
  struct endpt   *jointo[2];
  struct distrec *distbuf;

  used[0] = calloc(ncurves,sizeof(bool));
  used[1] = calloc(ncurves,sizeof(bool));

  jointo[0]  = calloc(ncurves,sizeof(struct endpt));
  jointo[1]  = calloc(ncurves,sizeof(struct endpt));

  distbuf = calloc(2*ncurves*2*ncurves,sizeof(struct distrec));
  
  int i,j,k,l,drcnt = 0;

  printf("joinvect: Building endpoint-endpoint distance array...");

  for(k=0;k<2;k++) {

    for(l=0;l<2;l++) {
      
      for(i=0;i<ncurves;i++) {
	
	for(j=i;j<ncurves;j++) {

	  if (j!=i || k != l) { /* We don't want to join a point to itself! */
	    
	    distbuf[drcnt].pts[0].curve = i;
	    distbuf[drcnt].pts[1].curve = j;
	    distbuf[drcnt].pts[0].sf = k;
	    distbuf[drcnt].pts[1].sf = l;
	    distbuf[drcnt].dist = torus_distance(
						k==0 ? nplc_start(curves[i],0) : nplc_end(curves[i],0),
						l==0 ? nplc_start(curves[j],0) : nplc_end(curves[j],0),
						twidth);
	  
	    drcnt++;

	  }

	}
	
      }

    }

  }

  printf("done.\n");

  printf("joinvect: Sorting %d point-distance records by distance...",drcnt);
  qsort(distbuf,drcnt,sizeof(struct distrec),distrec_cmp);
  printf("done.\n");

  printf("joinvect: Now joining points...");
  
  int bufcnt = 0;
  double maxdist = 0;
  
  for(i=0;i<ncurves;i++) { /* We must make ncurves links */

    for (;used[distbuf[bufcnt].pts[0].sf][distbuf[bufcnt].pts[0].curve] ||
	   used[distbuf[bufcnt].pts[1].sf][distbuf[bufcnt].pts[1].curve];bufcnt++); 

    jointo[distbuf[bufcnt].pts[0].sf][distbuf[bufcnt].pts[0].curve] = distbuf[bufcnt].pts[1];
    jointo[distbuf[bufcnt].pts[1].sf][distbuf[bufcnt].pts[1].curve] = distbuf[bufcnt].pts[0];

    used[distbuf[bufcnt].pts[0].sf][distbuf[bufcnt].pts[0].curve] = true;
    used[distbuf[bufcnt].pts[1].sf][distbuf[bufcnt].pts[1].curve] = true;

    if (distbuf[bufcnt].dist > maxdist) { maxdist = distbuf[bufcnt].dist; }

  }

  printf("done.\n");
  printf("joinvect: Maximum join distance is %g.\n",maxdist);

  /* We now reverse components as needed. */

  int  startc,thisc;
  bool reversed;
  bool *compused;

  compused = calloc(ncurves,sizeof(bool));

  for(startc=0;startc<ncurves;startc++) {

    if (!compused[startc]) {

      compused[startc] = true;
      order[startc] = jointo[1][startc].curve;
            
      int loopfail=0;

      for(thisc=jointo[1][startc].curve,reversed=(jointo[1][startc].sf==1);
	  thisc != startc;loopfail++) {
	
	compused[thisc] = true;
	order[thisc] = jointo[!reversed][thisc].curve;

	if (reversed) { 

	  nplc_reverse(curves[thisc]); 
	  
	}

	/* Now we increment for the next iteration of the loop. */

	if (reversed) {

	  reversed = (jointo[0][thisc].sf == 1);
	  thisc = jointo[0][thisc].curve;
	  
	} else {

	  reversed = (jointo[1][thisc].sf == 1);
	  thisc = jointo[1][thisc].curve;

	}

	if (loopfail == 10*ncurves) { 

	  printf("joinvect: Fatal error in order list generation.\n");
	  exit(1);

	}
	
      }

    }
    
  }
  
  free(compused);
  free(used[0]);
  free(used[1]);
  free(jointo[0]);
  free(jointo[1]);
  free(distbuf);

}

      
	
