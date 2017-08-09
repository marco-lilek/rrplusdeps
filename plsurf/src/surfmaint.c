/**********************************************/
/*            Surface Maintenance             */
/**********************************************/

/* This code is part of the surface library. */

#include "plsurf.h"
#include "plsurf_internal.h"

int new_surface(int verts,int faces,plc_color def_col,surface *surf)

  /* Procedure allocates space for a new surface.
     
     Since we don't know the number of vertices in each face,
     the procedure assumes that the answer is 0, and allocates
     all the face pointers to NULL.

     Procedure returns 1 if surface creation is successful,
                       0 if surface creation fails.

     */

{
  int i;

  /* First, we do some obvious sanity-checking on the input data. */

  if ((verts < 0) || (faces < 0)) {
  
    fprintf(stderr,"new_surface: Cannot create surface with %d verts," 
	    "%d faces.\n", verts, faces);

    return 0;

  }

  if (!goodcolor(def_col)) {
  
    fprintf(stderr,"new_surface: (%f,%f,%f,%f) is not a valid color.\n"
	    "             All numbers must be between 0 and 1.\n",
	    def_col.r,def_col.g,def_col.b,def_col.alpha);

    return 0;

  }

  if (surf == NULL) {

    fprintf(stderr,"new_surface: surf pointer is NULL (should point to allocated plsurf).\n");
    return 0;
    
  }
                 
  /* Now we are free to proceed. */

  surf->verts = verts;
  surf->faces = faces;
     
  surf->vert_buf = calloc(verts+1,sizeof(plc_vector));
  surf->vtf      = calloc(faces+1,sizeof(int));
  surf->face_buf = calloc(faces+1,sizeof(int*));
  surf->col_buf  = calloc(faces+1,sizeof(plc_color));

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

void new_face(surface *surf,int face_num,int vtf,...)

     /* Procedure allocates and creates a new face, given the face number,
	number of vertices in the face, and a list of vertex numbers which are
	specified using the va_arg protocol. 

        The missing arguments are expected to be int's indexing into the vertex
        buffer of surf. */
{
  va_list ap;
  int     i;
  
  /* We first check our input for sanity. */

  if (surf == NULL) {

    fprintf(stderr,"new_face: Can't create a face in a NULL surface!\n");
    exit(2);

  }

  if (!face_legal(surf,face_num)) {

    fprintf(stderr,"new_face: Face number %d is illegal in a surface with %d faces.\n",
	    face_num,surf->faces);

    exit(2);

  }

  if (vtf < 2) {

    fprintf(stderr,"new_face: Number of vertices %d is not permitted. \n",vtf);
    exit(2);

  }

  /* Now we can do the work. */

  surf->vtf[face_num] = vtf;
  surf->face_buf[face_num] = calloc(vtf,sizeof(int));

  va_start(ap,vtf);

  for(i=0;i<vtf;i++) {

    surf->face_buf[face_num][i] = va_arg(ap,int);

  }

  va_end(ap);
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

void reindex_face(surface *surf,int face,int start_vert)

     /* Procedure reorders the vertices in face, so that */
     /* start_vert comes first. Of course, this presumes */
     /* that start_vert _is_ a vertex of face! */

{
  int *scratch_buf,i,pos = {-1};

  for(i=0;i<surf->vtf[face];i++) {

    if (surf->face_buf[face][i] == start_vert) {

      pos = i;

    }

  }

  if (pos == -1) {

    /* We didn't find it! */

    fprintf(stderr,"reindex_face: Vertex %d is not on face %d.\n",
	    start_vert,face);

    return;

  }

  scratch_buf = calloc(surf->vtf[face],sizeof(int));

  for(i=0;i<surf->vtf[face];i++) {

    scratch_buf[i] = surf->face_buf[face][(i + pos) % surf->vtf[face]];

  }

  free(surf->face_buf[face]);
  surf->face_buf[face] = scratch_buf;

}

#ifdef converted

void triangulate_surface(surface *surf)

     /* Procedure triangulates the surface surf. */
     /* This procedure is destructive to surf, but */
     /* preserves as much information as possible. */
{
  int new_faces = 0;
  int i,j,this_face;
  int **triangles;
  plc_color NICE_BLUE = {0.1,0.1,0.8,1};
  polyline edge,planar;
  vector normal;

  int    *new_vtf;
  int   **new_face_buf;
  plc_color  *new_col_buf;

  /* We start by counting the new faces. */

  for(i=0;i<surf->faces;i++) {
  
    new_faces += surf->vtf[i] - 2;

  } 

  /* Now we make new buffers, with the required number of faces. */

  new_vtf      = calloc(new_faces,sizeof(int));
  new_face_buf = calloc(new_faces,sizeof(int*));
  new_col_buf  = calloc(new_faces,sizeof(plc_color));

  /* We are now ready to triangulate each face. */

  for(i=0,this_face=0;i<surf->faces;i++) {

    /* Debugging code. */

    /* fprintf(stderr,"face %d\n",i); */

    if(surf->vtf[i] > 3) {

      /* First, extract the outline of the face as a planar polyline... */

      face_outline(surf,i,&edge);
      normal = face_normal(surf,i);
      to_planar(&edge,&normal,&planar);

      /* Now triangulate it... */

      triangulate_planar_polyline(&planar,&(triangles));

      /* Now read the new faces into the face buffer. */
	      
      for(j=0;j<surf->vtf[i] - 2;j++) {

	new_vtf[this_face] = 3;
	new_face_buf[this_face] = calloc(3,sizeof(int));

	new_face_buf[this_face][0] = surf->face_buf[i][triangles[j][0]];
	new_face_buf[this_face][1] = surf->face_buf[i][triangles[j][1]];
	new_face_buf[this_face][2] = surf->face_buf[i][triangles[j][2]];

	new_col_buf[this_face] = surf->col_buf[i];

	this_face++;

      }

      /* Now clear the triangles buffer. */

      free(triangles);

    } else { 

      /* If the face is already a triangle... copy it.*/

      new_vtf[this_face] = 3;
      new_face_buf[this_face] = calloc(3,sizeof(int));

      new_face_buf[this_face][0] = surf->face_buf[i][0];
      new_face_buf[this_face][1] = surf->face_buf[i][1];
      new_face_buf[this_face][2] = surf->face_buf[i][2];
	
      new_col_buf[this_face] = surf->col_buf[i];

      this_face++;

    }

  }

  /* We now swap in the new data. */

  free(surf->face_buf);
  free(surf->vtf);
  free(surf->col_buf);

  surf->faces = new_faces;
  surf->vtf = new_vtf;
  surf->face_buf = new_face_buf;
  surf->col_buf = new_col_buf;

  /* We now need to kill all the lists. */

  kill_edge_list(surf);
  kill_incidence_lists(surf);

}

void subdivide_edges(surface *surf,int num_edges,int *edges)

     /* This procedure subdivides the given edges of surface.

	We expect num_edges to contain the number of edges to 
	operate on, and edges to be a pointer to an array of 
	these edge numbers.

	The procedure iterates through the list of edges, adding
	vertices to a new vert_buf, and rewriting the faces 
	in question as each edge is changed.

	This is all meaningless without the edge list, of course.
	We expect the edges list to contain no duplicates.

     */

{
  int i,j,k,l,m;
  vector *new_vert_buf;
  int new_verts;
  int *face_swap,face;
  int v[2];

  /* Step 0: Check to see that the edge list exists. */

  if (surf->edge_list == NULL || surf->edges == NOT_COMPUTED) {

    fprintf(stderr,"subdivide_edges: Edge list not computed.\n");
    exit(1);

  }

  /* Step 1: Allocate the new vertex buffer. */

  new_verts = surf->verts + num_edges;
  new_vert_buf = calloc(new_verts,sizeof(vector));

  memmove(new_vert_buf,surf->vert_buf,surf->verts*sizeof(vector));

  /* Step 2: Iterate through the list of edges. */

  for(j=0,i=surf->verts;j<num_edges;j++) {

    /* Debugging code. */

    if (edges[j] < 0 || edges[j] > surf->edges-1) {

      fprintf(stderr,"subdivide_edges: Illegal edge index (%d) at edge %d.\n",
	      edges[j],j);

    }

    /* First, write the new vertex into the new buffer. */

    v[0] = surf->edge_list[edges[j]][1];
    v[1] = surf->edge_list[edges[j]][2];

    midpoint(&surf->vert_buf[v[0]],
	     &surf->vert_buf[v[1]],
	     &new_vert_buf[i]);

    /* Now insert the vertex into each face on the list. */

    for(k=0;k<surf->edge_list[edges[j]][0];k++) {

      /* We make a new vertex list for this face. */

      face = surf->edge_list[edges[j]][2*k+3];

      /* Debugging code. */

      if (face < 0 || face > surf->faces-1) {

	fprintf(stderr,"subdivide_edges: Illegal face index %d.\n",
		face);

      }

      /* Now go on. */

      face_swap = calloc(surf->vtf[face]+1,sizeof(int));

      /* We search for the right combination in the face.  */
      /* The simplest way to do this is best. We reindex   */
      /* the face so that v[1] is first (if pos. oriented) */
      /* and v[0] is first (if neg. oriented). */

      if (surf->edge_list[edges[j]][2*k+4] == 1) {

	reindex_face(surf,face,v[1]);

      } else {

	reindex_face(surf,face,v[0]);

      }

      face_swap[0] = i;
	
      for(l=0;l<surf->vtf[face];l++) {

	face_swap[l+1] = surf->face_buf[face][l];

      }
	
      /* Now we swap the updated vertex list for the old one. */

      free(surf->face_buf[face]);
      surf->face_buf[face] = face_swap;
      surf->vtf[face]++;
      face_swap = NULL;

    }

    /* Now increment the vertex counter. */

    i++;

  }

  /* Step 3: Cleanup. */

  /* We have now subdivided all the edges in question,       */
  /* and updated the face vertex lists to reflect this fact. */
  /* At this point, we need to swap the new vertex buffer for the old. */

  free(surf->vert_buf);
  surf->vert_buf = new_vert_buf;
  surf->verts = new_verts;

  /* This destroys the incidence lists, so we kill them. */

  surf->edges = NOT_COMPUTED;
  free(surf->edge_list);

  kill_incidence_lists(surf);

}

void consol_recurser(surface *surf,int num_points,int *points,
		     int *lut,vector *ll,vector *ur,
		     double precision,int *this_newpoint)

     /* Procedure takes the vertices in the array points,
	which are to lie in the box given by ll and ur,
        and starts to look for duplicate vertices.

	If the size of the bbox is less than precision, 
	we stop the recursion, and identify these points
	in the lookup table lut. These are set to 
	this_newpoint, which is then incremented, and
	passed back up the chain.

	If not, we subdivide the box into 8 smaller cubes,
	binsort the points into eight smaller lists, and
	recurse to those cubes. 

        As a last act, we free the points list. */

{
  vector box_size;

  int    *sublist[2][2][2];
  int    sl_size[2][2][2] = {{{0,0},{0,0}},{{0,0},{0,0}}};

  int    i,j,k,l;
  vector new_ll,new_ur;

  /* First, we check the size of the box. */

  vectordiff(ur,ll,&box_size);
  
  if (box_size.x < precision && 
      box_size.y < precision &&
      box_size.z < precision) {

    /* We're done. We simply need to fill in 
       the appropriate values in the lookup table. */

    for(i=0;i<num_points;i++) {

      lut[points[i]] = (*this_newpoint);

    }

    (*this_newpoint)++;
    free(points);

    return;

  } 
  
  /* If the box is still too big, we subdivide. */
  /* The first step is allocating the little lists. */

  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      for(k=0;k<2;k++) {

	sublist[i][j][k] = calloc(num_points,sizeof(int));

      }
    }
  }

  /* Now we split up the points. */

  for(i=0;i<num_points;i++) {

    j = (surf->vert_buf[points[i]].x > ll->x + 0.5*box_size.x);
    k = (surf->vert_buf[points[i]].y > ll->y + 0.5*box_size.y);
    l = (surf->vert_buf[points[i]].z > ll->z + 0.5*box_size.z);

    sublist[j][k][l][sl_size[j][k][l]] = points[i];
    sl_size[j][k][l]++;

  }

  /* Now we recurse. */

  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      for(k=0;k<2;k++) {

	if (sl_size[i][j][k] > 0) {

	  new_ll.x = ll->x + i*0.5*box_size.x;
	  new_ll.y = ll->y + j*0.5*box_size.y;
	  new_ll.z = ll->z + k*0.5*box_size.z;

	  new_ur.x = new_ll.x + 0.5*box_size.x;
	  new_ur.y = new_ll.y + 0.5*box_size.y;
	  new_ur.z = new_ll.z + 0.5*box_size.z;

	  consol_recurser(surf,sl_size[i][j][k],
			  sublist[i][j][k],lut,
			  &new_ll,&new_ur,
			  precision,this_newpoint);
	}

      }
    }
  }

  /* Now we free the list handed to us and return. */

  free(points);

}

#endif 

void    subdivide_faces(surface *surf)

     /* Each face of surf must be a triangle. We perform a 
	barycentric subdivision of the surface, dividing each
	face into four smaller triangles. 

        We will triangulate the surface if we need to. */

{
  surface new_surf;
  int i,j,ns_face,os_face,*edges;

  /* First, we check to see that the edge and incidence lists are 
     computed for the surface we start with. */


  if (surf->edges == NOT_COMPUTED) {

    make_edge_list(surf);

  }

  if (surf->ftv == NULL) {

    make_incidence_lists(surf);

  }

  if (!is_triangulated(surf)) {

    fprintf(stderr,"subdivide_faces: Only works for triangulated surfaces.\n");
    exit(1);

  }

  /* Now we are ready to create the new surface. */

  new_surface(surf->verts + surf->edges,4*surf->faces,SURFACE_BLUE,&new_surf);

  /* First, we add the new vertices. */
  
  for(i=0;i<surf->verts;i++) {
    
    new_surf.vert_buf[i] = surf->vert_buf[i];

  }

  for(j=0;j<surf->edges;j++,i++) {

    new_surf.vert_buf[i] = plc_vlincomb(0.5,surf->vert_buf[surf->edge_list[j][1]],
					0.5,surf->vert_buf[surf->edge_list[j][2]]);

  }

  /* Now we work on the faces. It will be convenient to allocate */
  /* the face buffer and set vtf in one step at this point. */

  for(i=0;i<new_surf.faces;i++) {

    new_surf.vtf[i]      = 3;
    new_surf.face_buf[i] = calloc(3,sizeof(int));

  }

  for(os_face=0,ns_face=0;os_face<surf->faces;os_face++) {

    edges_on_face(surf,os_face,&edges);
    
    /* First, we fill in the three triangles around the perimeter. */

    for(i=0;i<3;i++) {
      
      new_surf.face_buf[ns_face][0] = surf->face_buf[os_face][i];
      new_surf.face_buf[ns_face][1] = edges[i] + surf->verts;
      new_surf.face_buf[ns_face][2] = edges[(i+2)%3] + surf->verts;

      ns_face++;

    }

    /* Now we fill in the central triangle. */

    new_surf.face_buf[ns_face][0] = edges[0] + surf->verts;
    new_surf.face_buf[ns_face][1] = edges[1] + surf->verts;
    new_surf.face_buf[ns_face][2] = edges[2] + surf->verts;

    ns_face++;

    /* Now, we free the edges array. */

    free(edges);

  }

  /* We have now completed the new surface. It remains to kill the old,
     and replace it with the new data. */

  kill_surface(surf);
  *surf = new_surf;

}

#ifdef CONVERTED

surface consolidate_surface(surface *surf,double precision)

     /* Procedure identifies all vertices in surf closer than
	<precision>, and produces a new surf, with fewer vertices,
	which is then returned. */

{
  surface result;
  int     i,j;
  int     *lut,*points;
  int     new_verts = {0};
  int     iterations;
  vector  ll,ur;

  /* First, we allocate the look-up table and the points table. */

  lut = calloc(surf->verts,sizeof(int));
  points = calloc(surf->verts,sizeof(int));

  for(i=0;i<surf->verts;i++) points[i] = i;
  
  /* Next, we identify groups of vertices in the old surface. */

  bbox_surface(surf,&ll,&ur);
  
  consol_recurser(surf,surf->verts,points,lut,
		  &ll,&ur,precision,&new_verts);

  /* Now we allocate the new surface. */

  new_surface(new_verts,surf->faces,SURFACE_BLUE,&result);

  /* Now we go ahead and fill it. */

  for(i=0;i<surf->verts;i++) {

    result.vert_buf[lut[i]] = surf->vert_buf[i];

  }

  for(i=0;i<surf->faces;i++) {

    result.vtf[i] = surf->vtf[i];
    result.face_buf[i] = calloc(surf->vtf[i],sizeof(int));
    
    for(j=0;j<surf->vtf[i];j++) {

      result.face_buf[i][j] = lut[surf->face_buf[i][j]];

    }

  }

  /* We are done with the lookup table and can free it. */

  free(lut);

  /* Last, we return the resulting surface. */

  return result;

}

void orient_recurser(surface *surf,int face,int *visited_faces)

     /* Procedure checks the surrounding faces for an edge-adjacent
	face previously visited by the algorithm. When it finds one,
	it alters the orientation of this face to match, then calls
	itself on all other neighboring faces. */

{
  int *faces,i,j,*edges,*temp_buf;

  /* First, we check the dumb stuff. */

  if (!face_legal(surf,face)) {

    fprintf(stderr,"orient_recurser: Illegal face number <%d>.\n",face);
    exit(1);

  }

  if (visited_faces[face] == TRUE) {

    return;

  }

  /* Now we can work! */

  adjacent_faces_to_face(surf,face,&faces);
  edges_on_face(surf,face,&edges);

  /* We are using a bit of forbidden knowledge here... we know
     that these buffers will match in that the face over edge 
     edges[i] should be faces[i]. */

  for(i=0;i<surf->vtf[face];i++) {

    if (faces[i] == -1) { 	/* If we are at a boundary, just continue. */

      continue;

    }

    if (visited_faces[faces[i]] == TRUE) {

      /* We now check to see if our orientations match. */

      if (surf->edge_list[edges[i]][4] ==
	  surf->edge_list[edges[i]][6]) {

	/* We must reverse the orientation of this face! */

	temp_buf = calloc(surf->vtf[face],sizeof(int));

	for(j=0;j<surf->vtf[face];j++) {

	  temp_buf[j] = surf->face_buf[face][surf->vtf[face] - j - 1];

	}

	for(j=0;j<surf->vtf[face];j++) {

	  surf->face_buf[face][j] = temp_buf[j];
	  
	}

	free(temp_buf);

	/* Now we must update all the edge records.         */
	/* We cannot simply recompute the edge list, since  */
	/* this will be too damn slow.                      */

	for(j=0;j<surf->vtf[face];j++) {

	  surf->edge_list[edges[j]]
	    [face_pos_in_edge_rec(surf,face,edges[j]) + 1] *= -1;
	  
	}

      }

    }

  }

  /* We have altered the orientation of this face (if we are going to!). */

  visited_faces[face] = TRUE;

  /* We now recurse to the adjacent faces. */

  for(i=0;i<surf->vtf[face];i++) {

    if (faces[i] != -1) { /* If we're at a boundary face */
                          /* the -1 will signal that. */

      orient_recurser(surf,faces[i],visited_faces);

    }

  }

  /* We now free the memory we've used here. */

  free(edges);
  free(faces);
  
}

void orient_surface(surface *surf)

     /* Procedure attempts to orient surf. This may not work.

	The procedure is smart enough to handle:

	non-orientable surfaces : We quite after surf->faces steps.
	disconnected surfaces   : We orient each orientable piece.

	The procedure is destructive to surf.

     */

{
  int i,*visited_faces;

  /* It's worth checking to see if we really HAVE to do this! */
  
  if (surf->edges == NOT_COMPUTED) {

    make_edge_list(surf);

  }

  if (is_oriented(surf)) {

    return;

  }
  
  /* Now we'll have to do some work. */
	
  visited_faces = calloc(surf->faces,sizeof(int));

  for(i=0;i<surf->faces;i++) {

    visited_faces[i] = FALSE;

  }

  /* Now we can start to orient. */

  for(i=0;i<surf->faces;i++) {

    if (!visited_faces[i]) {

      orient_recurser(surf,i,visited_faces);

    }

  }
  
  /* Last, we free our memory. */

  free(visited_faces);
  
  /* And reset the is_oriented flag on the surface! */

  surf->is_oriented = NOT_COMPUTED;
  
}

#endif
