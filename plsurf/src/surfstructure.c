/*********************************************/
/*     Structure of the Triangulation        */
/*********************************************/

/* This code is part of the surfaces library. */

#include"plsurf.h"
#include"plsurf_internal.h"


int edge_shared_faces(surface *surf,int edge_one,int edge_two,int **result)

     /* Result is set to a buffer containing the faces containing */
     /* both edges. The number of faces in result is returned. */

{
  int i,j,found;
  int *scratch;

  /* We simply search the incidence list for each edge, */
  /* noting our results in scratch as we find them. */

  /* We start by checking that the edge list is allocated. */

  if (surf->edge_list == NULL || surf->edges == 0) {

    fprintf(stderr,"edge_shared_faces: Edge list not allocated.\n");
    exit(1);

  }

  scratch = calloc(surf->faces,sizeof(int));
  found = 0;

  for(i=0;i<surf->edge_list[edge_one][0];i++) {

    for(j=0;j<surf->edge_list[edge_two][0];j++) {

      if (surf->edge_list[edge_one][2*i + 3] ==
	  surf->edge_list[edge_two][2*j + 3]) {

	scratch[found] = surf->edge_list[edge_one][2*i + 3];
	found++;

      }

    }

  }

  *result = calloc(found,sizeof(int));
  
  for(i=0;i<found;i++) {

    (*result)[i] = scratch[i];

  }

  free(scratch);

  return found;

}
	
int vert_shared_faces(surface *surf,int vert_one,int vert_two,int **result)

     /* Same as edge_shared faces, this procedure makes a list of faces */
     /* containing both verts, and puts it in result. The number of faces */
     /* in the list is returned. */

{
  int i,j;
  int *scratch, found;

  /* We simply search the incidence list for each vertex, */
  /* noting our results in scratch as we find them. */

  /* We start by checking that the incident_faces list is allocated. */

  if (surf->incident_faces == NULL) {

    fprintf(stderr,"vert_shared_faces: Incidence lists not allocated.\n");
    exit(1);

  }

  scratch = calloc(surf->faces,sizeof(int));
  found = 0;

  for(i=0;i<surf->ftv[vert_one];i++) {

    for(j=0;j<surf->ftv[vert_two];j++) {

      if (surf->incident_faces[vert_one][i] ==
	  surf->incident_faces[vert_two][j]) {

	scratch[found] = surf->incident_faces[vert_one][i];
	found++;

      }

    }

  }

  *result = calloc(found,sizeof(int));
  
  for(i=0;i<found;i++) {

    (*result)[i] = scratch[i];

  }

  free(scratch);

  return found;

}

int edge_vert_shared_faces(surface *surf,int edge,int vert,int **result)

     /* Same as above, this procedure makes a list of faces */
     /* containing both edge and vert, and puts it in result. */
     /* The number of faces in the list is returned. */

{
  int i,j;
  int *scratch, found;

  /* We simply search the incidence list for each vertex, */
  /* noting our results in scratch as we find them. */

  /* We start by checking that the incident_faces list is allocated. */

  if (surf->incident_faces == NULL || surf->edge_list == NULL) {

    fprintf(stderr,"edge_vert_shared_faces: Incidence or edge list "
	    "not allocated.\n");
    exit(1);

  }

  scratch = calloc(surf->faces,sizeof(int));
  found = 0;

  for(i=0;i<surf->edge_list[edge][0];i++) {

    for(j=0;j<surf->ftv[vert];j++) {

      if (surf->edge_list[edge][2*i + 3] ==
	  surf->incident_faces[vert][j]) {

	scratch[found] = surf->edge_list[edge][2*i + 3];
	found++;

      }

    }

  }

  *result = calloc(found,sizeof(int));
  
  for(i=0;i<found;i++) {

    (*result)[i] = scratch[i];

  }

  free(scratch);

  return found;

}

bool vert_on_face(surface *surf,int vert,int face)

     /* Procedure returns TRUE if vert is on face. */

{
  int j;
  bool found = false;

  for(j=0;j<surf->vtf[face];j++) {

    if (surf->face_buf[face][j] == vert) {

      found = true;

    }

  }

  return found;

}

bool edge_on_face(surface *surf,int edge,int face) 

     /* Procedure returns TRUE if edge is on face. */

{
  int j;
  bool found = false;

  /* We need the edge list for this. */

  if (surf == NULL) {

    fprintf(stderr,"edge_on_face: NULL surface passed.\n");
    exit(1);

  }

  if (surf->edge_list == NULL || surf->edges == 0) {

    fprintf(stderr,"edge_on_face: Edge list not computed.\n");
    exit(1);

  }

  for(j=0;j<surf->edge_list[edge][0];j++) {

    if(surf->edge_list[edge][2*j + 3] == face) {
    
      found = true;

    }

  }

  return found;

}

int vert_pos_on_face(surface *surf,int vert,int face) 

     /* Returns the position of vert on the face. */

{
  int i;

  /* First, we check the dumb stuff. */

  if (surf == NULL) {

    fprintf(stderr,"vert_pos_on_face: Surf pointer is null.\n");
    exit(1);

  }

  if (vert < 0 || vert >= surf->verts) {

    fprintf(stderr,"vert_pos_on_face: Vertex number %d is out of legal range"
	    " (0,%d) for surf.\n",vert,surf->verts);

    exit(1);

  }

  if (face < 0 || face >= surf->faces) {

    fprintf(stderr,"vert_pos_on_face: Face number %d is out of legal range"
	    " (0,%d) for surf.\n",face,surf->faces);

    exit(1);
    
  }

  /* Now we can actually work! */

  for(i=0;i<surf->vtf[face];i++) {

    if (surf->face_buf[face][i] == vert) {

      return i;

    }

  }

  fprintf(stderr,"vert_pos_on_face: Vert %d is not on face %d.\n",
	  vert,face);

  return -1;

}

int least_vert_pos_on_face(surface *surf,int face)

     /* Procedure searches the vertices on face, and returns the */
     /* position of the vertex with the least index. */

{
  int best_yet,best_pos,i;

  best_yet = surf->face_buf[face][0];
  best_pos = 0;

  for(i=0;i<surf->vtf[face];i++) {

    if (best_yet > surf->face_buf[face][i]) {

      best_yet = surf->face_buf[face][i];
      best_pos = i;

    }

  }

  return best_pos;

}

int  face_pos_at_vert(surface *surf,int face,int vert)

     /* Procedure returns the index of face in the incidence */
     /* list of faces meeting at vertex vert. If face is not */
     /* incident to vert, this is a fatal error. */

{
  int i;

  /* Step 0. Check everything. */

  if (!face_legal(surf,face) || !vert_legal(surf,vert)) {

    fprintf(stderr,"face_pos_at_vert: Illegal face or vertex index.\n");
    exit(1);

  }

  if (surf->incident_faces == NULL) { /* We need the incidence lists. */

    make_incidence_lists(surf);

  }

  /* Step 1. Search for face in list. */

  for(i=0;i<surf->ftv[vert];i++) {

    if (surf->incident_faces[vert][i] == face) {

      return i;

    }

  }

  fprintf(stderr,"face_pos_at_vert: Face %d is not incident to vert %d.\n",
	  face,vert);
  exit(1);

}     

int  face_pos_in_edge_rec(surface *surf,int face,int edge)

     /* Procedure searches the edge record for <edge>, looking
	for some mention of <face>. If it finds that <edge> is
	on <face>, we return the index in the surf->edge_list
	buffer corresponding to this face.

	If <edge> isn't on <face>, we return NOT_COMPUTED. */

{
  int i;

  /* We check the dumb stuff first. */

  if (surf->edges == NOT_COMPUTED) {

    fprintf(stderr,"face_pos_in_edge_list: Edge list not computed.\n");
    exit(1);

  }

  if (!face_legal(surf,face) || !edge_legal(surf,edge)) {

    fprintf(stderr,"face_pos_in_edge_list: face <%d> or edge <%d> illegal.\n",
	    face,edge);
    exit(1);

  }
  
  /* Now we can work. */

  for(i=0;i<surf->edge_list[edge][0];i++) {

    if (surf->edge_list[edge][2*i + 3] == face) {

      return 2*i + 3;

    }

  }

  return NOT_COMPUTED;

}

int  next_vert_on_face(surface *surf,int vert,int face)

     /* Procedure returns the next vertex in the list for this face,
        wrapping around to the beginning of the list if needed. */

{
  int vert_pos;

  /* As always, we check the dumb stuff. */

  if (!vert_on_face(surf,vert,face)) {
    
    fprintf(stderr,"next_vert_on_face: Vert <%d> is not on face <%d>.\n",vert,face);
    exit(1);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"next_vert_on_face: Face <%d> is out of range.\n",face);
    exit(1);

  }

  if (!vert_legal(surf,vert)) {

    fprintf(stderr,"next_vert_on_face: Vert <%d> is out of range.\n",vert);
    exit(1);

  }

  /* Now we can work! */

  vert_pos = vert_pos_on_face(surf,vert,face);
  return (vert_pos + 1) % surf->vtf[face];
  
}

int  vert_pos_opposite_edge_on_face(surface *surf,int edge,int face)

     /* Procedure returns the position of the vertex opposite edge
	on the (triangular) face given by face. */

{
  int pos[2];

  /* First, we check the dumb stuff. */

  if (!edge_legal(surf,edge) || !face_legal(surf,face)) {

    fprintf(stderr,"vert_pos_opposite_edge_on_face: Bad edge or face number "
	    "(%d,%d).\n",edge,face);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"vert_pos_opposite_edge_on_face: Edge %d is not on face %d.\n",
	    edge,face);

    exit(1);
    
  }

  /* We have now checked out the input data, and are ready to work. */

  pos[0] = vert_pos_on_face(surf,surf->edge_list[edge][1],face);
  pos[1] = vert_pos_on_face(surf,surf->edge_list[edge][2],face);

  /* We get the last formula from a trick of mod 3 arithmetic: Given
     two values in 0 1 2, we can get the third by computing

        2 * (0 + 1) = 2      mod 3
	2 * (0 + 2) = 4 = 1  mod 3
	2 * (1 + 2) = 6 = 0  mod 3
	
  */

  return ((2 * (pos[0] + pos[1])) % 3);

}

void face_edges_meeting_at_vert(surface *surf,int face,int vert,int *edges)

     /* Procedure returns the edge numbers of the two edges meeting
	at vert on face, in the order given by the orientation of the
	face. 

	If vert is not on face, we return an error and quit. 

        We expect edges to point to a buffer of at least two integers. */
{
  int pos;
  int vtf;
  
  /* First, we have to check the dumb stuff. */

  if (surf == NULL) {

    fprintf(stderr,"face_edges_meeting_at_vert: NULL surface passed.\n");
    exit(1);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"face_edges_meeting_at_vert: Face number %d illegal.\n",
	    face);
    exit(1);

  }

  if (!vert_legal(surf,vert)) {

    fprintf(stderr,"face_edges_meeting_at_vert: Vert number %d illegal.\n",
	    vert);
    exit(1);

  }

  if (!vert_on_face(surf,vert,face)) {

    fprintf(stderr,"face_edges_meeting_at_vert: Vert %d isn't on face %d.\n",
	    vert,face);
    exit(1);

  }

  /* Now we can get to work! */

  pos = vert_pos_on_face(surf,vert,face);
  vtf = surf->vtf[face];

  edges[0] = edge_number(surf,surf->face_buf[face][(pos + vtf - 1) % vtf],
			 surf->face_buf[face][pos]);

  edges[1] = edge_number(surf,surf->face_buf[face][pos],
			 surf->face_buf[face][(pos + 1) % vtf]);

}


bool faces_adjacent(surface *surf,int face_one,int face_two,int *edge)
     
  /* This procedure checks to see whether these faces share an edge. 
     If so, the edge number is returned. */

{
  int i,j;
  bool found_one = false,found_two = false,adj = false;

  /* We search the edge list. This is slow, but we're going to have */
  /* to do it sometime, so we might as well do it now. */

  for(i=0;i < surf->edges;i++) {

    found_one = found_two = false;

    for(j=0;j < surf->edge_list[i][0];j++) {

      if (surf->edge_list[i][2*j + 3] == face_one) {

	found_one = true;

      }

      if (surf->edge_list[i][2*j + 3] == face_two) {

	found_two = true;

      }

    }

    if (found_one && found_two) {

      adj = true;
      *edge = i;

    }

  }

  return adj;

}

static int edge_comp(const void *a,const void *b)

     /* Procedure compares two edge_list entries, returning */
     /* -1 if a < b, 0 if a == b, and 1 if a > b. */

{
  int *A,*B;

  if (a == NULL || b == NULL) {

    fprintf(stderr,"edge_comp: Passed NULL pointer to compare.\n");
    exit(2);

  }

  A = *((int **)(a));
  B = *((int **)(b));

  if (A[1] < B[1]) {

    return -1;

  } else if (A[1] > B[1]) {

    return 1;

  } else if (A[2] < B[2]) {

    return -1;

  } else if (A[2] > B[2]) {

    return 1;

  } else {

    return 0;

  }
}

int  faces_vertex_adjacent(surface *surf,int f_one,int f_two,int **verts)

     /* Procedure decides whether f_one and f_two are vertex-adjacent. 
	If so, returns the number of shared vertices, and sets *verts to
	a pointer to the number of vertices in question. */

{
  int i,j,found;
  int *shared;

  /* First, we check the indices. */

  if (!face_legal(surf,f_one) || !face_legal(surf,f_two)) {

    fprintf(stderr,"faces_vertex_adjacent: Face index out of range.\n");
    exit(1);

  }

  /* Now we work. First, we allocate a big buffer. */

  shared = calloc(surf->vtf[f_one] + surf->vtf[f_two],sizeof(int));
  found  = 0;

  /* Now we search for shared vertices. */

  for(i=0;i<surf->vtf[f_one];i++) {

    for(j=0;j<surf->vtf[f_two];j++) {

      if (surf->face_buf[f_one][i] == surf->face_buf[f_two][j]) {

	shared[found++] = surf->face_buf[f_one][i];

      }

    }

  }

  if (found > 0) {		/* If we found anything, we return it. */

    (*verts) = shared;
    return found;

  } else {			/* If not, we free buffer & quit. */

    free(shared);
    return 0;

  }

}

void adjacent_faces_to_face(surface *surf,int face,int **faces) 

     /* Procedure allocates and returns a list of the edge-adjacent
	faces to face. If, for some odd reason, a face shares two 
	edges with <face>, it will be listed twice in this list. */

{
  int i,*edges;
  
  /* First, we check the dumb stuff. */

  if (surf->edges == NOT_COMPUTED) {

    make_edge_list(surf);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"adjacent_faces_to_face: Face <%d> is out of range.\n",
	    face);
    exit(1);

  }

  /* Now we can work in peace. First, we get a list of the edges on face.*/

  edges_on_face(surf,face,&edges);

  /* Now we make a buffer */

  *faces = calloc(surf->vtf[face],sizeof(int));
  
  /* Last, we look over each edge to the adjacent face! */

  for(i=0;i<surf->vtf[face];i++) {

    (*faces)[i] = adjacent_face(surf,edges[i],face);

  }

  /* Now we can free the edge list. */

  free(edges);

}


int edge_number(surface *surf,int vert_one,int vert_two)

     /* Procedure searches the edge list for the given edge.     */
     /* If not found, returns -1 as the edge number.             */
     /* Since the edge list is kept in dictionary order, we      */
     /* use bsearch to find the number. */

{
  int v[3];
  int **result;
  int *vptr;

  /* Our first step is to put the edge into normal form. */
  /* To do so, we make a "fake" edge_list entry... */

  v[0] = 0;

  if (vert_one == vert_two) {

    return -1;

  } else if (vert_one < vert_two) {

    v[1] = vert_one; v[2] = vert_two;

  } else {

    v[1] = vert_two; v[2] = vert_one;

  }

  /* Now we can search... */

  vptr = &(v[0]);

  result = bsearch(&vptr,&(surf->edge_list[0]),surf->edges,
		   sizeof(int *),&edge_comp);

  if (result == NULL) {

    return -1;

  } else {

    return (int)(result - surf->edge_list);

  }

}

int adjacent_face(surface *surf,int edge,int face)

     /* Procedure returns the face number of the face adjacent to */
     /* face over edge. If edge is not on face, or if there is no */
     /* other face over edge, or if more than one face is adjacent, */
     /* an error is reported, and -1 is returned. */

{

  /* First, we check for the really bad case- when surf is not 
     a top manifold. */

  if (surf->edge_list[edge][0] > 2) {

    fprintf(stderr,"adjacent_face: Edge %d is on %d faces.\n",
	    edge,surf->edge_list[edge][0]);

    return -1;

  } 

  /* Next, we check that surf->edge_list[edge][0] is at least one.
     If not, the data structure has been corrupted, so we print an
     error message, and terminate the program. */

  if (surf->edge_list[edge][0] < 1) {

    fprintf(stderr,"adjacent_face: Illegal edge record at edge %d.\n",edge);
    exit(1);

  }

  /* Last, we consider the case where edge is a boundary edge. 
     If so, we silently return -1 to indicate that there is no 
     adjacent edge. */

  if (surf->edge_list[edge][0] == 1) {

    return -1;

  }

  if (surf->edge_list[edge][3] == face) {

    return surf->edge_list[edge][5];

  }

  if (surf->edge_list[edge][5] == face) {

    return surf->edge_list[edge][3];

  }

  fprintf(stderr,"adjacent_face: Edge %d not on face %d.\n",
	  edge,face);

  return -1;

}

void edge_positions_on_face(surface *surf,int edge,int face,int edgepos[2])

     /* Procedure computes the positions of the vertices in edge on face,
	and returns them, IN THE ORDER GIVEN BY THE ORIENTATION OF THE FACE. 
	Thus, the first entry in edgepos may correspond to the second vertex
	of the face, if the edge is negatively oriented on face. */

{
  int swap;

  /* First, we check the dumb stuff. */

  if (!edge_legal(surf,edge) || !face_legal(surf,face)) {

    fprintf(stderr,"edge_positions_on_face: Bad edge or face number (%d,%d).\n",
	    edge,face);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"edge_positions_on_face: Edge %d is not on face %d.\n",
	    edge,face);

    exit(1);
    
  }

  /* We have now checked out the input data, and are ready to work. */

  edgepos[0] = vert_pos_on_face(surf,surf->edge_list[edge][1],face);
  edgepos[1] = vert_pos_on_face(surf,surf->edge_list[edge][2],face);

  /* We want to arrange these so that they are consecutive (mod 3). */
  /* Either this is true, and 

                   edgepos[1] - 1 = edgepos[0] mod 3, 

     or this is false, and

                   edgepos[0] - 1 = edgepos[1] mod 3.

     We test for the second case. */

  if (((edgepos[0] - 1 + 3) % 3) == edgepos[1]) {

    swap = edgepos[0];
    edgepos[0] = edgepos[1];
    edgepos[1] = swap;

  }

}
  
int  edge_orientation_on_face(surface *surf,int edge,int face)

     /* Procedure returns +1/-1 to indicate the orientation of
	EDGE on FACE. If EDGE is not on FACE, we print an error
	and quit. */

{
  int fpos;

  /* As always, first we check the dumb stuff. */

  if (surf == NULL) {

    fprintf(stderr,"edge_orientation_on_face: NULL surface passed.\n");
    exit(1);

  }

  if (!edge_legal(surf,edge)) {

    fprintf(stderr,"edge_orientation_on_face: Edge number %d not legal.\n",
	    edge);
    exit(1);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"edge_orientation_on_face: Face number %d not legal.\n",
	    face);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"edge_orientation_on_face: Edge %d not on face %d.\n",
	    edge,face);
    exit(1);

  }

  /* Now we're ready to work. */

  fpos = face_pos_in_elist_record(surf,face,edge);

  return surf->edge_list[edge][fpos+1];

}
	
void edges_on_face(surface *surf,int face,int **edges)

     /* Procedure makes a list of the edges on a particular face.   */
     /* If, for some perverse reason, an edge occurs twice, we list */
     /* it twice in the result. */

{
  int i;

  /* First, we check the dumb stuff. */

  if (surf->edges == NOT_COMPUTED) {

    make_edge_list(surf);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"edges_on_face: Face <%d> is out of range.\n",face);
    exit(1);

  }

  /* Now we can work in peace. */
  
  *edges = calloc(surf->vtf[face],sizeof(int));
  
  for(i=0;i<surf->vtf[face];i++) {
    
    (*edges)[i] = edge_number(surf,
			      surf->face_buf[face][i],
			      surf->face_buf[face][(i+1)%surf->vtf[face]]);

  }
}

int  next_edge_on_face(surface *surf,int face,int edge)

     /* Procedure works counterclockwise around the verts of the given
	face, listing the edge next to the given edge IN THE ORDER OF
	THE ORIENTATION OF THE FACE. */

{
  int edgepos[2];
  int vert_one,vert_two;
  
  /* First, we check the stupid stuff. */

  if (surf->edges == NOT_COMPUTED) {

    make_edge_list(surf);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"next_edge_on_face: Face <%d> is out of range.\n",face);
    exit(1);

  }

  if (!edge_legal(surf,edge)) {

    fprintf(stderr,"next_edge_on_face: Edge <%d> is out of range.\n",edge);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"next_edge_on_face: Edge <%d> is not on face <%d>.\n",edge,face);
    exit(1);

  }

  /* Now we are ready to get to work. */

  edge_positions_on_face(surf,edge,face,edgepos);
  vert_one = edgepos[1];
  vert_two = next_vert_on_face(surf,vert_one,face);

  return edge_number(surf,vert_one,vert_two);

}
  
void edges_at_vert(surface *surf,int vert,int *edges)

     /* Procedure computes the edges coming in to a vertex,
	and writes them, in counterclockwise order, to the buffer
	edges, which is assumed to be allocated and at least
	as large as surf->ftv[vert]. 

        Since ftv is incremented by one when the vertex is a 
        boundary vert, we always return surf->ftv[vert] edges,
        assuming that surf is a topological manifold. */

{
  int i,face,pos;

  /* First, we check the dumb stuff. */

  if (!vert_legal(surf,vert)) {
    
    fprintf(stderr,"edges_at_vert: Vertex %d illegal.\n",vert);
    exit(1);

  }

  if (surf->edges == NOT_COMPUTED) {

    make_edge_list(surf);

  }

  if (surf->ftv == NULL) {

    make_incidence_lists(surf);

  }

  /* Now we can get to work. */

  for(i=0;i<surf->ftv[vert];i++) {

    face = surf->incident_faces[vert][i];
   
    if (face != NO_FACE) {   /* We're not at the end of boundary vert. */

      pos  = vert_pos_on_face(surf,vert,face);
      edges[i] = edge_number(surf,
			     surf->face_buf[face][pos],
			     surf->face_buf[face][(pos+1)%surf->vtf[face]]);

    } else {		     /* We need to take care of the last edge.  */

      face = surf->incident_faces[vert][surf->ftv[vert]-2];
      pos  = vert_pos_on_face(surf,vert,face);

      edges[i] = edge_number(surf,
			     surf->face_buf[face][pos],
			     surf->face_buf[face][(pos-1+surf->vtf[face])%surf->vtf[face]]);

    }

  }
			   
}


