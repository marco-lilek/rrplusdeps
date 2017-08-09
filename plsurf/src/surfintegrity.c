/**********************************************/
/*         Surface Integrity Checks           */
/**********************************************/

/* This code is part of the surface library.  */

#include"plsurf.h"
#include"plsurf_internal.h"

#ifdef CONVERTED

bool sanity_surface(surface *surf)

     /* Procedure runs some basic checks on the surface. */
     /* Returns TRUE silently if surface passes.         */
     /* Returns false with an error if surface fails.    */

{
  int surfgood = {true};

  if (!edge_list_ok(surf)) {

    fprintf(stderr,"sanity_surface: Edge list is not ok.\n");
    surfgood = false;

  }

  if (!faces_ok(surf)) {

    fprintf(stderr,"sanity_surface: Face list is not ok.\n");
    surfgood = false;

  }

  if (!colors_good(surf)) {

    fprintf(stderr,"sanity_surface: Color list is not good.\n");
    surfgood = false;

  }

  return surfgood;

}

#endif

bool is_triangulated(surface *surf)

     /* Procedure checks whether surf is triangulated. */

{
  int face;
  bool flag = true;

  for(face = 0;face < surf->faces;face++) {

    if (surf->vtf[face] != 3) {

      flag = false;

    }

  }

  return flag;

}

#ifdef CONVERTED

bool is_top_mfld(surface *surf)

  /* Procedure checks whether surf represents a */
  /* topological (pl) manifold with boundary.   */

  /* Eventually, this will actually work. At the */
  /* moment, we only implement the first step of */
  /* the algorithm, checking that each edge occurs */
  /* in only 1 or 2 faces. */

{
  int i;

  if (surf->is_top_mfld != NOT_COMPUTED) {

    return surf->is_top_mfld;

  } else {

    /* We'll need the edge list to make this work. */

    if (surf->edges <= 0 || surf->edge_list == NULL) {

      make_edge_list(surf);

    }

    /* Now we check the number of faces each edge appears in. */

    for(i=0;i<surf->edges;i++) {

      if (!((surf->edge_list[i][0] == 1) || (surf->edge_list[i][0] == 2))) {
      
	fprintf(stderr,"is_top_mfld: Edge %d appears in %d faces. \n"
	             "             Surface is not a topological manifold.\n\n",
		i,surf->edge_list[i][0]);

	surf->is_top_mfld = false;
	return false;

      }

    }

    /*  fprintf(stderr,"is_top_mfld: Warning! Surface has passed first-stage "
	"tests only!\n\n"); */

    surf->is_top_mfld = true;
    return true;
    
  }

}

bool is_oriented(surface *surf)

  /* If surface is a topological manifold, we can ask whether it carries */
  /* a consistent orientation at the moment. We do so by inspecting the  */
  /* edge list, asking that each edge occuring in two faces appears with */
  /* opposite orientations. */

     /* We return true or false silently, as this is used in practice */
     /* before applications of "orient". */

{
  int i;

  if (surf->is_oriented != NOT_COMPUTED) {

    return surf->is_oriented;

  } else {

    if (!is_top_mfld(surf)) {

      fprintf(stderr,"is_oriented: Warning! Cannot evaluate orientation \n"
	      "             of a non-manifold. Returning false.\n");
     
      return false;
      
    }

    for (i=0;i<surf->edges;i++) {

      if (surf->edge_list[i][0] == 2) {

	if (surf->edge_list[i][4] + surf->edge_list[i][6] != 0) {

	  surf->is_oriented = false;
	  return false;

	}

      }

    }

    surf->is_oriented = true;

    return true;

  }

}

bool is_closed(surface *surf)

  /* Procedure checks whether the surface */
  /* is closed (without boundary).        */

  /* We note that the boundary is defined */
  /* whether or not surf is a topological */
  /* manifold, so we do not check for this. */

{
  int i;

  if (surf->is_closed != NOT_COMPUTED) {

    return surf->is_closed;

  } else {

    /* First, we make sure that the edge list has been computed. */

    if (surf->edges <= 0 || surf->edge_list == NULL) {

      make_edge_list(surf);

    }

    /* Now we check the # of faces each edge occurs in. */

    for(i=0;i<surf->edges;i++) {

      if (surf->edge_list[i][0] < 2) {

	surf->is_closed = false;
	return false;

      }

    }

    surf->is_closed = true;
    return true;
    
  }

}

bool edge_list_ok(surface *surf) 

     /* Procedure checks the edge list to make sure that all references 
	are legal, and that the edge list looks ok. */

{
  int j,k;

  if (surf->edges == NOT_COMPUTED) {

    return true;

  }

  if (surf->edges <= 0) {

    return false;

  }

  for(j=0;j<surf->edges;j++) {

    if (surf->edge_list[j][0] <= 0) {

      fprintf(stderr,"edge_list_ok: Bad number of faces %d at edge %d.\n",
	      surf->edge_list[j][0],j);

      return false;

    }

    if (surf->edge_list[j][1] < 0 ||
	surf->edge_list[j][1] > surf->verts-1) {

      fprintf(stderr,"edge_list_ok: Bad vertex number %d at edge %d.\n",
	      surf->edge_list[j][1],j);

    }

    if (surf->edge_list[j][2] < 0 ||
	surf->edge_list[j][2] > surf->verts-1) {

      fprintf(stderr,"edge_list_ok: Bad vertex number %d at edge %d.\n",
	      surf->edge_list[j][2],j);

    }

    for(k=0;k<surf->edge_list[j][0];k++) {

      if (surf->edge_list[j][2*k+3] < 0 ||
	  surf->edge_list[j][2*k+3] > surf->faces-1) {

	fprintf(stderr,"edge_list_ok: Bad face number %d at edge %d.\n",
		surf->edge_list[j][2*k+3],j);

      }

      if (!(surf->edge_list[j][2*k+4] == 1 ||
	    surf->edge_list[j][2*k+4] == -1)) {

	fprintf(stderr,"edge_list_ok: Bad orientation %d at edge %d.\n",
		surf->edge_list[j][2*k+3],j);

      }  

    }

  }

  return true;

}

bool verts_unique(surface *surf)

  /* Procedure checks the vertex list to make sure that it does not */
  /* contain any repeating vertices. */

{
  int i,j;

  for(i=0;i<surf->verts-1;i++) {

    for(j=i+1;j<surf->verts;j++) {

      if (segmentlength(&surf->vert_buf[i],&surf->vert_buf[j]) < EPSILON) {

	fprintf(stderr,"verts_unique: Vertex %d (%f,%f,%f) and \n"
		       "              vertex %d (%f,%f,%f) \n"
		       "              of surface appear to be duplicates.\n\n",
		i,surf->vert_buf[i].x,surf->vert_buf[i].y,surf->vert_buf[i].z,
		j,surf->vert_buf[j].x,surf->vert_buf[j].y,surf->vert_buf[j].z);

	return false;

      }

    }

  }
  
  return true;

} 

#endif

bool faces_ok(surface *surf)

  /* Procedure checks whether the faces of surf are all ok. */

{
  int i;

  for(i=0;i<surf->faces;i++) {

    if (!face_ok(surf,i)) {

      fprintf(stderr,"faces_ok: Warning! Face %d of surface is not ok.\n\n",
	      i);

      return false;

    }

  }

  return true;

}

#ifdef CONVERTED

bool faces_planar(surface *surf) 

  /* Procedure checks the faces of surf for planarity. */

{
  int i;

  for(i=0;i<surf->faces;i++) {

    if (!face_planar(surf,i)) {

      fprintf(stderr,"faces_planar: WARNING! Face %d of surface "
	      "isn't planar.\n\n",i);

      return false;

    }

  }

  return true;

}

bool colors_good(surface *surf)

  /* Procedure checks that each of the colors in surf's col_buf */
  /* refer to legal GEOMVIEW rgba colors. */

{
  int i;

  for (i=0;i<surf->faces;i++) {

    if (!goodcolor(surf->col_buf[i])) {

      fprintf(stderr,"colors_good: WARNING! Color of face %d is bad!\n"
	             "             (%f,%f,%f,%f).\n\n",
	      i,
	      surf->col_buf[i].r,
	      surf->col_buf[i].g,
	      surf->col_buf[i].b,
	      surf->col_buf[i].a);

      return false;

    }

  }

  return true;

}
	    
#endif 

bool face_ok(surface *surf,int i) 

  /* Procedure decides whether the ith face of surf consists of */
  /* nonrepeating vertex references in the range 0..verts-1. */

{
  int j,k;

  for(j=0;j<surf->vtf[i]-1;j++) {

    for(k=j+1;k<surf->vtf[i];k++) {

      if ((surf->face_buf[i][j] < 0) || 
	  (surf->face_buf[i][j] > surf->verts-1)) {

	fprintf(stderr,"face_ok: Vertex reference %d of face %d of surface \n"
		"is out of the legal range 0..%d.\n\n",j,i,surf->verts-1);
		
	return false;

      }

      if (surf->face_buf[i][j] == surf->face_buf[i][k]) {

	fprintf(stderr,"face_ok: Vertex reference %d occurs in position %d "
		"and %d of face %d.\n\n",surf->face_buf[i][j],j,k,i);

	return false;

      }

    }

  }

  return true;

}

#ifdef CONVERTED

bool face_planar(surface *surf,int i)

  /* Procedure decides whether the ith face of surf lies in a single plane. */

{
  int j;
  vector normal;
  vector *diffs;

  /* First, if the face is a triangle, we're done. */

  if (surf->vtf[i] == 3) {
    
    return true;

  }

  /* If not, we must compute. Our method is to compute */
  /* a surface normal by moving the origin to the first */
  /* vertex, then taking cross products. */

  diffs = calloc(surf->vtf[i],sizeof(vector));

  for(j=0;j<surf->vtf[i];j++) {

    vectordiff(&surf->vert_buf[surf->face_buf[i][j]],
	       &surf->vert_buf[surf->face_buf[i][0]],
	       &diffs[j]);
  }

  cross(&diffs[1],&diffs[2],&normal);

  /* Now we check that this normal vector is perpendicular to all vertices. */

  for(j=3;j<surf->vtf[i];j++) {

    if (fabs(dot(&diffs[j],&normal)) > 
	(0.0001)*norm(&diffs[j])*norm(&normal)) {

      return false;

    }

  }

  return true;

}

#endif

bool goodcolor(plc_color col)

  /* Procedure checks the values of a color structure. */

{
  if (col.r < 0 || col.r > 1) {

    return false;

  }

  if (col.g < 0 || col.g > 1) {

    return false;

  }

  if (col.b < 0 || col.b > 1) {

    return false;

  }

  if (col.alpha < 0 || col.alpha > 1) {

    return false;

  }

  return true;

}

bool vert_legal(surface *surf,int vert)

     /* Procedure checks to see whether vert is within bounds for surf. */

{

  if (vert < 0 || vert > surf->verts-1) {

    return false;

  }

  return true;

}

bool edge_legal(surface *surf,int edge) 
     
     /* Procedure checks to see whether edge is within bounds for surf. */
     /* Of course, before we do so, we must check that the edge list has */
     /* been computed for surf. */

{

  if (surf->edges == NOT_COMPUTED || surf->edge_list == NULL) {

    return false;
    
  }

  if (edge < 0 || edge > surf->edges-1) {

    return false;

  }

  return true;

}

bool face_legal(surface *surf,int face) 

     /* Procedure checks to see whether face is within bounds for surf. */

{

  if (face < 0 || face > surf->faces - 1) {

    return false;

  } 

  return true;

}
