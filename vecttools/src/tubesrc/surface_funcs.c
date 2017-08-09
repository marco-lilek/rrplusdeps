/*

    surface_funcs: Contains utility functions for the modified surface type. 

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "plCurve.h"

#include"tube.h"
#include"maxmin.h"

int new_surface(int verts,int faces,surface *surf)

  /* Procedure allocates space for a new surface.
      
     Since we don't know the number of vertices in each face,
     the procedure assumes that the answer is 0, and allocates
     all the face pointers to NULL.

     Procedure returns 1 if surface creation is successful,
                       0 if surface creation fails.

     */

{
  int i;
  plc_color def_col = {1.0,1.0,1.0,1.0};

  assert(surf != NULL);
  assert(faces > 0);
  assert(verts > 0);

  /* Now we are free to proceed. */

  surf->verts = verts;
  surf->faces = faces;
     
  surf->vert_buf = (plc_vector *)(calloc(verts+1,sizeof(plc_vector)));
  surf->vtf      = (int *)(calloc(faces+1,sizeof(int)));
  surf->face_buf = (int **)(calloc(faces+1,sizeof(int*)));
  surf->col_buf  = (plc_color *)(calloc(faces+1,sizeof(plc_color)));

  if ((surf->vert_buf == NULL) || (surf->vtf == NULL) || 
      (surf->face_buf == NULL) || (surf->col_buf == NULL)) {

    fprintf(stderr,"new_surface: Not enough memory to allocate surface"
	    " with %d vertices and %d faces.\n",verts,faces);

    return 0;

  }

  /* We now set sensible default values for the remaining entries. */

  for (i=0;i<faces;i++) {

    surf->col_buf[i] = def_col;

  }

  surf->edges = NOT_COMPUTED;
  surf->edge_list = NULL;
  surf->incident_faces = NULL;
  surf->incident_angles = NULL;
  surf->ftv = NULL;

  return 1;
}


void write_surf_to_OFF(surface *surf,FILE *outfile)

  /* Procedure writes a surface to a Geomview-compatible OFF file. */
  /* We do use the keyword OFF at the beginning of the file, and   */
  /* we always provide 4 digit real color information. */

{
  int i,j;

  fprintf(outfile,"OFF \n");
  fprintf(outfile,"%d %d %d\n",surf->verts,surf->faces,surf->edges);
  
  for(i=0;i<surf->verts;i++) {

    fprintf(outfile,"%f %f %f \n",
	    surf->vert_buf[i].c[0],
	    surf->vert_buf[i].c[1],
	    surf->vert_buf[i].c[2]);

  }

  for(i=0;i<surf->faces;i++) {

    fprintf(outfile,"%d ",surf->vtf[i]);
    
    for(j=0;j<surf->vtf[i];j++) {

      fprintf(outfile,"%d ",surf->face_buf[i][j]);

    }

    fprintf(outfile,"%f %f %f %f \n",
	    surf->col_buf[i].r,
	    surf->col_buf[i].g,
	    surf->col_buf[i].b,
	    surf->col_buf[i].alpha);

  }

}


void kill_surface(surface *surf)

  /* Procedure frees the memory associated with a surface. */

{
  int i;

  /* First, we free all the individual face vertex lists. */

  for(i=0;i<surf->faces;i++) {

    free(surf->face_buf[i]);
  
  }

  /* Now, if the edge_list has been filled, we free that. */

  if (surf->edge_list != NULL && surf->edges > 0) {

    for(i=0;i<surf->edges;i++) {

      free(surf->edge_list[i]);

    }

  }
  
  /* Now, we kill the incidence lists, if they've been filled. */

  kill_incidence_lists(surf); 

  /* Next, we free the vtf list, the face, vertex, and color buffers. */

  free(surf->vtf);      surf->vtf = NULL;
  free(surf->face_buf); surf->face_buf = NULL;
  free(surf->vert_buf); surf->vert_buf = NULL;
  free(surf->col_buf);  surf->col_buf = NULL;
  free(surf->edge_list); surf->edge_list = NULL;

  /* Last, we reset the verts, faces, edges, and area. */

  surf->verts = 0;
  surf->faces = 0;
  surf->edges = NOT_COMPUTED;

}

void kill_edge_list(surface *surf) 

     /* Procedure kills the edge list. */

{
  int j;

  if (surf->edges > 0) {

    for(j=0;j<surf->edges;j++) {

      free(surf->edge_list[j]);

    }

    free(surf->edge_list);
    surf->edges = NOT_COMPUTED;

  }
}

void kill_incidence_lists(surface *surf)

/* Procedure kills the incident_faces and incident_angles lists. */

{
  int i;

  if (surf->incident_faces != NULL && surf->verts != 0) {

    for(i=0;i<surf->verts;i++) {

      free(surf->incident_faces[i]);

    }

    free(surf->incident_faces);
    surf->incident_faces = NULL;

  }  

  if (surf->incident_angles != NULL && surf->verts != 0) {

    for(i=0;i<surf->verts;i++) {

      free(surf->incident_angles[i]);

    }

    free(surf->incident_angles);
    surf->incident_angles = NULL;

  }  

  if (surf->ftv != NULL) {

    free(surf->ftv);
    surf->ftv = NULL;

  }

}


void build_ring(int numsteps, double radius, 
		plCurve *L,
		int cp,
		int vt,
		//plc_vector loc, 
		plc_vector frameA, 
		plc_vector frameB,
		surface *tube,
		double s,
		struct uvbuf *uvb)

/* Procedure creates a list of vertices in the tube structure for a ring 
   of the tube. We assume that tube has space allocated for these verts
   and that tube->verts is set to the index of the first free vertex 
   in the buffer. 

   We set (s,theta) coordinates for the vertices on the ring as well.

   If the turning angle at this vertex is more than 5 degrees, we actually
   generate an elliptical ring which takes into account that this ring is 
   not perpendicular to the incident edges. */

{
  plc_vector loc;
  double theta;
  int i;
  double twopi = 4*acos(0);
  double turnangle,fivedeg = (2*acos(0))*(5.0/180.0);
  double norangle = 0,norrad = radius;
  bool   ok;
  static int ringnum = 0;
  plc_vector tan,nor,bin;
  
  loc = L->cp[cp].vt[vt];
  
  turnangle = plc_turning_angle(L,cp,vt,&ok);

  if (!(turnangle < fivedeg)) {

      /* We will need to generate an elliptical cross-section
	 to make the tube look right. We first compute the normal
	 vector of the curve and write its projection on the plane
	 in terms of the frame vectors. */

    plc_vector ftan,rtan;
    bool ok;
    
    tan = tube_tangent(L,cp,vt,&ok);
    assert(ok);
    ftan = plc_normalize_vect(
			      plc_vect_diff(L->cp[cp].vt[vt+1],
					    L->cp[cp].vt[vt]), &ok);

    rtan = plc_normalize_vect(
			      plc_vect_diff(L->cp[cp].vt[vt-1],
					    L->cp[cp].vt[vt]), &ok);
    nor 
      = plc_normalize_vect(plc_vect_sum(ftan,rtan),&ok);
    /* Set nor to the vector in the ftan, tan plane normal to tan. */

    bin = plc_normalize_vect(plc_cross_prod(nor,tan),&ok);

    assert(fabs(plc_dot_prod(nor,ftan) - plc_dot_prod(nor,rtan)) < 1e-12);

    /* We now assert that nor should be in the frameA, frameB plane */

    assert(fabs(pow(plc_dot_prod(nor,frameA),2.0) + pow(plc_dot_prod(nor,frameB),2.0) - 1.0) < 1e-12);
    assert(fabs(plc_norm(frameA) - 1) < 1e-12);
    assert(fabs(plc_norm(frameB) - 1) < 1e-12);

    /* Now we continue */

    norrad = radius/plc_norm(plc_cross_prod(nor,ftan));
    norangle = atan2(norrad*plc_dot_prod(frameB,nor),radius*plc_dot_prod(frameA,nor));

  }
  
  for(i=0,theta=0;i<numsteps;theta+=(twopi)/numsteps,i++) {
    
    if (turnangle < fivedeg) {
      
      tube->vert_buf[tube->verts] = 
	plc_vlincomb(radius*cos(theta),frameA,
		     radius*sin(theta),frameB);
    
    } else {

      double angle;

      angle = -(theta-norangle); /* The effective angle on the ellipse */

      tube->vert_buf[tube->verts] = 
	plc_vlincomb(norrad*cos(angle),nor,
		     radius*sin(angle),bin);
      
    }
	
    plc_M_add_vect(tube->vert_buf[tube->verts],loc);

    if (MAKERINGS) {

      rings->cp[ringnum].vt[i] = tube->vert_buf[tube->verts];

    }

    if (uvb != NULL) {

      uvb->uv[tube->verts].c[1] = theta;
      uvb->uv[tube->verts].c[0] = s;

    }

    (tube->verts)++;
  }

  ringnum++;

}
  

void match_rings(int startA, int endA, int startB, int endB, 
		 int ab_ofs, plc_color ringcolor, surface *tube) 

     /* This procedure matches two rings of vertices on the surface <tube>
	using Bresenham's algorithm if the number of vertices is different. 

	We assume that the surface has space allocated for these triangles, 
	and that the number of faces in the surface data structure is set to
	the number of currently filled slots in the data structure. 

        We expect that the "A band" is before the "B band" on the core curve
        of the tube. We also expect that the regions of vertices are circular.
        The 0th vertex of the "A band" will match with the "ab_ofs" vertex of
        the "B band". */

{

  int thisface;
  int xbandstart,ybandstart,a,b,x0,y0;
  int e;
  int xofs, yofs;
  int vbA, vbB;

  /* First, we check the input for sanity. */

  assert(tube != NULL);
  assert(startA >= 0 && startA < endA && endA <= tube->verts);
  assert(startB >= 0 && startB < endB && endB <= tube->verts);
  
  /* Now we start to work. */

  vbA = endA - startA; vbB = endB - startB;

  /* We need to decide which of these intervals has more triangles. */
  /* By assumption, the "x-step" steps in the direction of the one 
     with more triangles */

  if (vbA > 1 || vbB > 1) {

    /* If there are no triangles here at all, we do nothing. */
    
    x0 = y0 = 0;
    thisface = tube->faces;

    if (vbA > vbB) {	    /* We are setting things so "line" is ax + by = 0 */
                            /* with line through (0,0) and (vbA,vbB) */
      xbandstart = startA;
      ybandstart = startB;
      
      a = vbB; 
      b = -vbA;

      xofs = 0;
      yofs = ab_ofs;
      
    } else { 			/* Same deal, except (0,0) and (vbB,vbA) */
      
      xbandstart = startB;
      ybandstart = startA;
      
      a = vbA;
      b = -vbB;

      xofs = ab_ofs;
      yofs = 0;
      
    }
    
    for(;x0<-b-1;) {
      
      e = 2*a*(x0+1) + 2*b*(y0) - 1;
      
      /* We are now (finally) ready to add faces. */
      /* We always take an x step. */
      
      tube->vtf[thisface] = 3;
      tube->face_buf[thisface] = (int *)(calloc(3,sizeof(int)));
      tube->col_buf[thisface] = ringcolor;
      
      if (xbandstart == startA) { /* Preserve a consistent orientation on the surface. */
	
	tube->face_buf[thisface][0] = xbandstart + ( (x0 + xofs) % vbA );
	tube->face_buf[thisface][1] = xbandstart + ( (x0 + xofs + 1) % vbA);
	tube->face_buf[thisface][2] = ybandstart + ( (y0 + yofs) % vbB);
	
      } else {
	
	tube->face_buf[thisface][0] = xbandstart + ( (x0 + xofs + 1) % vbB);
	tube->face_buf[thisface][1] = xbandstart + ( (x0 + xofs) % vbB);
	tube->face_buf[thisface][2] = ybandstart + ( (y0 + yofs) % vbA);
	
      }
      
      thisface++; x0++;
      
      if (e > 0 && y0 < a-1) { 		/* We take a y step as well   */
	                                /* Notice that doesn't happen */
	                                /* if y0 is already maxed */
	tube->vtf[thisface] = 3;
	tube->face_buf[thisface] = (int *)(calloc(3,sizeof(int)));
	tube->col_buf[thisface] = ringcolor;

	if (ybandstart == startA) {
	  
	  tube->face_buf[thisface][0] = ybandstart + ( (y0 + yofs) % vbA);
	  tube->face_buf[thisface][1] = ybandstart + ( (y0 + yofs + 1) % vbA);
	  tube->face_buf[thisface][2] = xbandstart + ( (x0 + xofs) % vbB);
	  
	} else {
	  
	  tube->face_buf[thisface][0] = ybandstart + ( (y0 + yofs + 1) % vbB);
	  tube->face_buf[thisface][1] = ybandstart + ( (y0 + yofs) % vbB);
	  tube->face_buf[thisface][2] = xbandstart + ( (x0 + xofs) % vbA);
	  
	}
	
	thisface++; y0++;
	
      }
      
    }    

    /* We now need to sew together the first and last vertices of the tubes. */
    
    tube->vtf[thisface] = 3;
    tube->face_buf[thisface] = (int *)(calloc(3,sizeof(int)));
    tube->col_buf[thisface] = ringcolor;
    
    if (xbandstart == startA) { /* Preserve a consistent orientation on the surface. */
      
      tube->face_buf[thisface][0] = xbandstart + ( ( x0 + xofs ) % vbA);
      tube->face_buf[thisface][1] = xbandstart + ( ( 0  + xofs ) % vbA);
      tube->face_buf[thisface][2] = ybandstart + ( ( 0  + yofs ) % vbB);
      
      /* Modified this to fix sawtooth bug. */

    } else {
      
      tube->face_buf[thisface][0] = xbandstart + ( ( 0 + xofs) % vbB);
      tube->face_buf[thisface][1] = xbandstart + ( (x0 + xofs) % vbB);
      tube->face_buf[thisface][2] = ybandstart + ( (y0 + yofs) % vbA);
      
    }
    
    thisface++;
    
    tube->vtf[thisface] = 3;
    tube->face_buf[thisface] = (int *)(calloc(3,sizeof(int)));
    tube->col_buf[thisface] = ringcolor;
    
    if (xbandstart == startA) { /* Preserve a consistent orientation */
      /* on the surface. */
      
      tube->face_buf[thisface][0] = ybandstart + ( ( 0 + yofs) % vbB);
      tube->face_buf[thisface][1] = ybandstart + ( (y0 + yofs) % vbB);
      tube->face_buf[thisface][2] = xbandstart + ( (x0 + xofs) % vbA);
      
    } else {
      
      tube->face_buf[thisface][0] = ybandstart + ( (y0 + yofs) % vbA);
      tube->face_buf[thisface][1] = ybandstart + ( (0  + yofs) % vbA);
      tube->face_buf[thisface][2] = xbandstart + ( (0  + xofs) % vbB);
      
    }
    
    thisface++;
    
    /* We have now finished creating all the triangles. */
    
    tube->faces = thisface;

  }

}
      

int compute_numsteps(plCurve *link, int cp, int vt, 
		     double rad, double stepratio)

/* This procedure computes the angular step used to generate 
   the tube mesh at this vertex. It is designed so that the 
   triangles in the mesh will be approximately square. There
   are many caveats in such a design-- one is that the minimum
   number of sides in a tube is set to 15.*/

{
  double lenIn, lenOut, corestep;
  plc_vector diff;
  int numsteps;					
  double twopi = 4*acos(0);

  diff = plc_vect_diff(link->cp[cp].vt[vt+1],link->cp[cp].vt[vt]);
  lenIn = plc_M_norm(diff);

  diff = plc_vect_diff(link->cp[cp].vt[vt],link->cp[cp].vt[vt-1]);
  lenOut = plc_M_norm(diff);

  if (lenIn < 1e-6) { corestep = lenOut; } 
  else if (lenOut < 1e-6) { corestep = lenIn; }
  else { corestep = (lenIn + lenOut)/2.0; }

  numsteps = intmax(2,(int)(ceil(stepratio*twopi*(rad/corestep))),15);

  if (numsteps < MINSIDES) { numsteps = MINSIDES; }

  return numsteps;
}

void close_frame(plCurve *L, int cmp, plCurve *frameA, plCurve *frameB, 
		 int numsteps, int *ofs)

/* Procedure alters the frame to make frameA->cp[cmp].vt[frameA->cp[cmp].nv-1] line up 
   with one of the numsteps vertices in the frameA->cp[cmp].vt[0], frameB->cp[cmp].vt[0]
   plane. It then returns the number of this vertex in "ofs". */

{
  double holonomy;
  double theta_step;
  double pi = 2.0*acos(0);
  double ip;

  plc_vector lastA,nfA,nfB;
  int    i;

  assert(L != NULL); assert(frameA != NULL); assert(frameB != NULL); assert(ofs != NULL);

  /* First, we compute the holonomy of the frame. */

  lastA = frameA->cp[cmp].vt[frameA->cp[cmp].nv-1]; /* This saves typing. */

  assert(fabs(plc_M_norm(frameA->cp[cmp].vt[0]) - 1) < 1e-11);
  assert(fabs(plc_M_norm(frameB->cp[cmp].vt[0]) - 1) < 1e-11);
  assert(fabs(plc_M_norm(lastA) - 1) < 1e-11);

  holonomy = atan2(plc_M_dot(lastA,frameB->cp[cmp].vt[0]),plc_M_dot(lastA,frameA->cp[cmp].vt[0]));

  if (holonomy < 0) holonomy += 2*pi;  

  /* Now we adjust the frame to make the holonomy a multiple of
     (2pi)/numsteps. */

  modf(holonomy/(2*pi/(double)(numsteps)), &ip);
  *ofs  = (int)(ip);

  theta_step  = holonomy - ip*(2*pi/(double)(numsteps)); 
  theta_step /= -(double)(L->cp[cmp].nv);

  for(i=0;i<L->cp[cmp].nv;i++) {


    plc_M_vlincomb(nfA,
		   cos(i*theta_step),frameA->cp[cmp].vt[i],
		   sin(i*theta_step),frameB->cp[cmp].vt[i]);

    /* octrope_vlincombine(frameA->cp[cmp].vt[i],cos(i*theta_step),
       frameB->cp[cmp].vt[i],sin(i*theta_step),
       nfA); */

    plc_M_vlincomb(nfB,
		   -sin(i*theta_step),frameA->cp[cmp].vt[i],
		   cos(i*theta_step),frameB->cp[cmp].vt[i]);

    /* octrope_vlincombine(frameA->cp[cmp].vt[i],-sin(i*theta_step),
       frameB->cp[cmp].vt[i],cos(i*theta_step),
       nfB); */

    frameA->cp[cmp].vt[i] = nfA;
    frameB->cp[cmp].vt[i] = nfB;

  }

}


surface *make_tube(plCurve *L, double (*radius)(int comp,int vert),
		   double stepratio, 
		   plCurve *frameA, plCurve *frameB,
		   int capped,
		   struct uvbuf *uvb)

{
  int total_verts;
  int total_faces, ringsteps[2], ringstart[2];
  int i,j;
  int ofs;

  int first_ringsteps = {0},first_ringstart = {0};
  plc_color ringcolor = {0,0,0,1};

  surface *tube_surf;

  /* We start by checking our pointers. */

  assert(L != NULL);  assert(frameA != NULL); assert(frameB != NULL);
  
  plc_fix_wrap(L); 
  plc_fix_wrap(frameA); 
  plc_fix_wrap(frameB);

  /* First, we compute the number of vertices and faces in the tube. */

  int *nv = NULL;
  int *cc = NULL;
  bool *open = NULL;
  int cp = 0;

  if (MAKERINGS) {

    nv = calloc(plc_num_verts(L),sizeof(int));
    cc = calloc(plc_num_verts(L),sizeof(int));
    open = calloc(plc_num_verts(L),sizeof(int));

  }

  for(i=0,total_verts=0;i<L->nc;i++) {

    for(j=0;j<L->cp[i].nv;j++) {

      int steps;

      steps = compute_numsteps(L,i,j,radius(i,j),stepratio);
      total_verts += steps;

      if (MAKERINGS) {
	
	open[cp] = false;
	nv[cp++] = steps;
	
      }

    }

  }

  if (MAKERINGS) {

    rings = plc_new(cp,nv,open,cc);

  }

  if (uvb != NULL) {

    uvb->verts = total_verts;
    uvb->uv = calloc(uvb->verts,sizeof(plc_vector));
    assert(uvb->uv != NULL);

  }
    
  /* We claim: total_faces <= 2 * total_verts. 

     Proof. Each ring of the tube contains one segment per vert.
     And each of these segments is the base of at most two triangular faces.
     None of these triangles is shared between rings, since such a 
     face would have at least 4 edges. 

     In fact, the total number of faces is less when the tube is 
    .open, since those ring edges only bound one face. */

  total_faces = 2 * total_verts;
  
  tube_surf = (surface *)(malloc(sizeof(surface)));

  new_surface(total_verts,total_faces + 2*L->nc,tube_surf);
  tube_surf->faces = 0;	   /* These are set to 0 to indicate to build_ring */
  tube_surf->verts = 0;    /* and match_rings the number of faces constructed */
                           /* so far. */

  /* We now create the tube itself. */

  assert(L->nc == frameA->nc);
  assert(L->nc == frameB->nc);

  double s;

  for(i=0;i<L->nc;i++) {

    assert((L->cp[i].nv == frameA->cp[i].nv) &&
	   (frameA->cp[i].nv == frameB->cp[i].nv));
    assert((L->cp[i].open == frameA->cp[i].open) && 
	   (frameA->cp[i].open == frameB->cp[i].open));  


    ringstart[0] = tube_surf->verts;
    ringsteps[0] = compute_numsteps(L,i,0,radius(i,0),stepratio);

    /* Now if this component is closed, we're going to have to modify
       the frame to close exactly, which it won't generally do. We
       also save the vertex numbers of the first ring, since we'll
       come back to it later. */

    if (!L->cp[i].open) {

      close_frame(L,i,frameA,frameB,ringsteps[0],&ofs);
      first_ringsteps = ringsteps[0];
      first_ringstart = ringstart[0];

    }
    
    s = 0;

    build_ring(ringsteps[0],radius(i,0),L,i,0,
	       frameA->cp[i].vt[0],frameB->cp[i].vt[0],
	       tube_surf,s,uvb);

    /* We generate a cap face for this tube if needed. */

    if (capped && L->cp[i].open) { 

      tube_surf->vtf[tube_surf->faces] = ringsteps[0];
      tube_surf->face_buf[tube_surf->faces] = 
	(int *)(malloc(ringsteps[0]*sizeof(int)));

      for(j=0;j<ringsteps[0];j++) {

	tube_surf->face_buf[tube_surf->faces][j] = 
	  ringstart[0] + j;

      }

      tube_surf->faces++;

    }

    for(j=1;j<L->cp[i].nv;
	j++,
	  ringsteps[0]=ringsteps[1],
	  ringstart[0]=ringstart[1]) {

      /* First, we set the color for this ring. */

      if (L->cp[i].cc == 0) {

	ringcolor.r = ringcolor.g = ringcolor.b = ringcolor.alpha = 1.0;

      } else if (L->cp[i].cc == 1) {

	ringcolor = L->cp[i].clr[0];

      } else if (L->cp[i].cc == L->cp[i].nv) {

	ringcolor = L->cp[i].clr[j];

      } else {

	assert(L->cp[i].cc == 0 || L->cp[i].cc == 1 || \
	       L->cp[i].cc == L->cp[i].nv);

	exit(1);

      }

      /* Now we build the ring. */

      ringsteps[1] = compute_numsteps(L,i,j,radius(i,j),stepratio);
      ringstart[1] = tube_surf->verts;

      assert(ringstart[1] + ringsteps[1] <= total_verts);

      s += plc_distance(L->cp[i].vt[j-1],L->cp[i].vt[j]);

      build_ring(ringsteps[1],radius(i,j),
		 L,i,j,frameA->cp[i].vt[j],frameB->cp[i].vt[j],
		 tube_surf,s,uvb);

      assert(ringsteps[0] + ringsteps[1] + tube_surf->faces <= total_faces);

      match_rings(ringstart[0],ringstart[0]+ringsteps[0],
		  ringstart[1],ringstart[1]+ringsteps[1],
		  0,
		  ringcolor,
		  tube_surf);

    }

    if (!L->cp[i].open) { /* The component is closed. */

      match_rings(ringstart[0],ringstart[0]+ringsteps[0],
		  first_ringstart,first_ringstart + first_ringsteps,
		  ofs,
		  ringcolor,
		  tube_surf);

    } else if (capped) { /* Generate a cap face. */

        tube_surf->vtf[tube_surf->faces] = ringsteps[0];
	tube_surf->face_buf[tube_surf->faces] = 
	  (int *)(malloc(ringsteps[0]*sizeof(int)));

	for(j=0;j<ringsteps[0];j++) {
	  
	  tube_surf->face_buf[tube_surf->faces][j] = 
	    ringstart[0] + ringsteps[0] - 1 - j;
	  
	}
	
	tube_surf->faces++;
	
    }
    
  }

  return tube_surf;
  
}

