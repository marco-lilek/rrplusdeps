/***********************************************/
/*      Edge and Incidence List Code           */
/***********************************************/

/*  This code is part of the surface library.  */

#include "plsurf.h"
#include "plsurf_internal.h"

struct edge_struct {

  int v[2];
  int orientation;
  int face;

};

struct edgelist_struct {

  int nfaces;
  int v[2];
  int *faces;
  int *orientations;

};

int edge_compare(const void *a,const void *b)

{
  const struct edge_struct *esa = (const struct edge_struct *) a;
  const struct edge_struct *esb = (const struct edge_struct *) b;
     
  if (esa->v[0] != esb->v[0]) {

    return esa->v[0] - esb->v[0];

  } else {

    return esa->v[1] - esb->v[1];

  }

}

void make_edge_list(surface *surf)

     /* Procedure compiles a list of the edges of the surface, and the */
     /* faces they appear in, with the following format:               */
  
     /* <#faces> <v_1> <v_2> <f_1> <o_1> <f_2> <o_2> ...               */
  
     /* where the v_i are the endpoints of the edge, the f_i are the   */
     /* faces the edge appears in, and the o_i, all equal to +1 or -1  */
     /* decide whether the edge appears as written, or reversed.       */

     /* The edge list should be in dictionary order when complete! */

{
  struct edge_struct *esbuf;

  /* Step 1. Compile a list of edges. */
  
  int esbuf_size,nfound=0;

  esbuf_size = 1000;
  esbuf = calloc(esbuf_size,sizeof(struct edge_struct));

  int i,j;
  struct edge_struct *esp;

  for(i=0;i<surf->faces;i++) {

    for(j=0;j<surf->vtf[i];j++) {

      esp = &esbuf[nfound];

      if (surf->face_buf[i][j] < surf->face_buf[i][(j+1)%surf->vtf[i]]) {
	
	esp->v[0] = surf->face_buf[i][j];
	esp->v[1] = surf->face_buf[i][(j+1)%surf->vtf[i]];
	esp->face = i;
	esp->orientation = 1;

      } else {

	esp->v[1] = surf->face_buf[i][j];
	esp->v[0] = surf->face_buf[i][(j+1)%surf->vtf[i]];
	esp->face = i;
	esp->orientation = -1;

      }

      nfound++;
      
      if (nfound == esbuf_size) {

	esbuf_size *= 2;
	esbuf = realloc(esbuf,esbuf_size*sizeof(struct edge_struct));

      }

    }

  }

  /* Step 2. Sort the list. */

  qsort(esbuf,nfound,sizeof(struct edge_struct),edge_compare);
  
  /* Step 3. Collapse duplicates and write the data. */

  surf->edge_list = calloc(nfound,sizeof(int *));
  surf->edges = 0;

  int nfaces;
  
  for(i=0;i<nfound;) {

    /* Count duplicate edges. */
    
    for(j=i;!edge_compare(&esbuf[i],&esbuf[j]);j++);
    nfaces = j-i;

    /* Now allocate space and fill in data. */

    surf->edge_list[surf->edges] = calloc(3 + 2*nfaces,sizeof(int));

    surf->edge_list[surf->edges][0] = nfaces;
    surf->edge_list[surf->edges][1] = esbuf[i].v[0];
    surf->edge_list[surf->edges][2] = esbuf[i].v[1];

    for(j=i;!edge_compare(&esbuf[i],&esbuf[j]);j++) {

      surf->edge_list[surf->edges][3+2*(j-i)] = esbuf[j].face;
      surf->edge_list[surf->edges][4+2*(j-i)] = esbuf[j].orientation;

    }

    /* Finally, increment surf->edges and i. */

    surf->edges++;
    i = j;
    
  }
    
}

int edges_surface(surface *surf)

  /* Procedure computes the total number of distinct edges in the surface */
  /* Procedure returns the value, and sets it in surf's data structure. */
  
{
  int i;
  int first_guess = 0;
  int this_face,this_vert,forw_face,forw_vert;
  int edge_dups,edge_first,edge_last,forw_first,forw_last;

  /* First, we check our assumption that the faces of surf are all ok. */

  if (!faces_ok(surf)) {

    fprintf(stderr,"edges_surface: Warning! Can't count edges, since not all "
	    "faces are ok.\n");

    return 0;

  }

  /* Next, we compute our first guess. */

  for(i=0;i<surf->faces;i++) {

    first_guess += surf->vtf[i];

  }

  /* Now, we isolate a particular edge. */

  for(this_face=0;this_face < surf->faces-1;this_face++) {

    for(this_vert=0;this_vert < surf->vtf[this_face];this_vert++) {

      /* We want to find the first and last vertices. */

      edge_first = surf->face_buf[this_face][this_vert];

      /* The last one is special if we are at the end of the list. */

      if (this_vert == surf->vtf[this_face]-1) {

	edge_last = surf->face_buf[this_face][0];

      } else {

	edge_last = surf->face_buf[this_face][this_vert + 1];

      }

      edge_dups = 0;

      /* We now want to compare edge_first, edge_last with other edges. */

      for(forw_face = this_face+1;forw_face < surf->faces;forw_face++) {

	for(forw_vert = 0;forw_vert < surf->vtf[forw_face];forw_vert++) {

	  forw_first = surf->face_buf[forw_face][forw_vert];

	  if (forw_vert == surf->vtf[forw_face]-1) {

	    forw_last = surf->face_buf[forw_face][0];

	  } else {

	    forw_last = surf->face_buf[forw_face][forw_vert+1];

	  }

	  /* We have now isolated a comparison edge. We check to */
	  /* see if it is the same as our original edge, remembering */
	  /* that it may have the opposite orientation. */

	  if (edge_first == forw_first && edge_last == forw_last) {

	    edge_dups++;

	  }

	  if (edge_first == forw_last && edge_last == forw_first) {

	    edge_dups++;

	  }

	}

      }

      /* We have now counted the number of times this edge shows up. */
      /* Since we will detect recounts again later, we simply subtract */
      /* one from our guess if edge_dups > 0. */

      if (edge_dups > 0) {
	
	first_guess--;

      }

    }

  }

  /* We're done! We set the edges field in surf, then return our findings. */

  surf->edges = first_guess;

  return first_guess;

}


void make_incidence_lists(surface *surf)
     
     /* Procedure makes a list of the faces incident to each vertex of */
     /* surf, in the format : */

     /* <f_1> ... <f_n> */

     /* and fills in the incident_faces structure in *surf. The procedure */
     /* also fills in the "incident_angles" data structure. The faces appear
        in counterclockwise order around the vertex. */

     /* This procedure requires the edge list to be computed. Hence,
	we check for the edge list, and if it is not present, compute
	it silently. */

     /* If the vertex is a boundary vertex, special rules apply. In these
	cases, we add an extra face f_{n+1}, with the value NO_FACE (= -1)
	to indicate that the record ends. 

	In these cases, incident_angles also receives an extra entry 
	(set to 0), and ftv is incremented as well. 

	These changes are picked up later by appropriate code elsewhere. */

{
  int this_vert,this_face,vert_idx,vtf;
  int included_faces,i;
  int ofs,shared_edge,new_face,prev_face;
  int *iscratch;
  double *ascratch;

  /* We check the basics first. */

  if (surf == NULL) {

    fprintf(stderr,"make_incidence_lists: surf pointer is NULL.\n");
    exit(1);

  }

  if (surf->edges == NOT_COMPUTED || surf->edge_list == NULL) {

    make_edge_list(surf);

  }

  /* We start by killing the previous versions of these lists. */

  kill_incidence_lists(surf);

  /* Now we make new ones. */ 

  surf->incident_faces = calloc(surf->verts,sizeof(int *));      // These are POINTERS
  surf->incident_angles = calloc(surf->verts,sizeof(double *));  // This, too.
  surf->ftv = calloc(surf->verts,sizeof(int));                   // But this is really an int.

  for(this_face=0;this_face < surf->faces;this_face++) {

    /* Our method is to iterate through the list of faces, proceeding
       around each vertex that we find, as we find it. If a vertex has
       already been handled, we'll ignore it and proceed to the next. */

     for(vert_idx=0;vert_idx < surf->vtf[this_face];vert_idx++) {

       this_vert = surf->face_buf[this_face][vert_idx];

       if (surf->incident_faces[this_vert] == NULL) {
	
	/* Since we don't know how many faces we'll find, we make some
	   scratch buffers to hold the results. */

	 iscratch = calloc(surf->faces,sizeof(int));
	 ascratch = calloc(surf->faces,sizeof(double));

	 /* We now initialize the scratch buffers with this face's data. */

	 included_faces = 0;
	 iscratch[included_faces] = this_face;
	 
	 ascratch[included_faces] = face_angle_at_vert(surf,
						       this_face,
						       this_vert);

	 /* We are now ready to move around to the next face. We do this
	    by finding the inbound edge to our vert, and locating the 
	    face adjacent over that edge. */

	 ofs = vert_pos_on_face(surf,this_vert,this_face);
	 vtf = surf->vtf[this_face];

	 shared_edge = edge_number(surf,
				   surf->face_buf[this_face]
				     [(ofs + vtf - 1)%vtf],
				   surf->face_buf[this_face][ofs]);

	 new_face = adjacent_face(surf,shared_edge,this_face);
	 
	 /* We also compute the previous face. */

	 shared_edge = edge_number(surf,
 				   surf->face_buf[this_face][ofs],
				   surf->face_buf[this_face][(ofs+1)%vtf]);

	 prev_face = adjacent_face(surf,shared_edge,this_face);

	 /* We now iterate through all such faces. */

	 for(included_faces = 1;
	     new_face != -1 && new_face != this_face;
	     included_faces++) {

	    iscratch[included_faces] = new_face;
	 
	    ascratch[included_faces] = face_angle_at_vert(surf,
							  new_face,
							  this_vert);

	    ofs = vert_pos_on_face(surf,this_vert,new_face);
	    vtf = surf->vtf[new_face];

	    shared_edge = edge_number(surf,
				      surf->face_buf[new_face]
				        [(ofs + vtf - 1)%vtf],
				      surf->face_buf[new_face][ofs]);

	    new_face = adjacent_face(surf,shared_edge,new_face);
	 
	 } 

	 /* We are now done. There are two cases here. */

	 if (new_face == -1 && prev_face == -1) {

	   /* We are at a boundary vertex, and we've covered the entire
	      span of incident faces. We need to add the extra record to
	      the end of incident_faces, and incident_angles, and copy
	      the lot to their final buffers. */

	   surf->incident_faces[this_vert] = calloc(included_faces + 1,
						    sizeof(int));

	   surf->incident_angles[this_vert] = calloc(included_faces + 1,
						     sizeof(double));

	   surf->ftv[this_vert] = included_faces + 1;

	   for(i=0;i<included_faces;i++) {

	     surf->incident_faces[this_vert][i] =  iscratch[i];
	     surf->incident_angles[this_vert][i] = ascratch[i];

	   }

	   surf->incident_faces[this_vert][i] = NO_FACE;
	   surf->incident_angles[this_vert][i] = 0.0;

	 }

	 if (new_face != -1) {

	   /* We are at an interior vertex, and have wrapped around. */

	   surf->incident_faces[this_vert] = calloc(included_faces,
						    sizeof(int));

	   surf->incident_angles[this_vert] = calloc(included_faces,
						     sizeof(double));

	   surf->ftv[this_vert] = included_faces;

	   for(i=0;i<surf->ftv[this_vert];i++) {

	     surf->incident_faces[this_vert][i] =  iscratch[i];
	     surf->incident_angles[this_vert][i] = ascratch[i];

	   }

	 }

	 free(iscratch);
	 free(ascratch);

      }

     }

  }

}

int  face_pos_in_elist_record(surface *surf,int face,int edge)

     /* Procedre returns the index of a face in the record 
	corresponding to an edge. If EDGE is not on FACE, we
	print an error and quit. */

{
  int i;
  
  /* As always, first we check the dumb stuff. */

  if (surf == NULL) {

    fprintf(stderr,"face_pos_in_elist_record: NULL surface passed.\n");
    exit(1);

  }

  if (!edge_legal(surf,edge)) {

    fprintf(stderr,"face_pos_in_elist_record: Edge number %d not legal.\n",
	    edge);
    exit(1);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"face_pos_in_elist_record: Face number %d not legal.\n",
	    face);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"face_pos_in_elist_record: Edge %d not on face %d.\n",
	    edge,face);
    exit(1);

  }

  /* Now we're ready to work. */

  for(i=3;i<surf->edge_list[edge][0]*2 + 2;i+=2) {

    if (surf->edge_list[edge][i] == face) {

      return i;

    }

  }

  /* Since we have already checked that edge was on face, we should
     never get here. If we do, there is a bug somewhere in the above
     code, so we print an error and terminate. */

  fprintf(stderr,"face_pos_in_elist_record: Internal error.\n");
  exit(1);

}


