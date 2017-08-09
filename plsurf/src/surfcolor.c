/************************************************/
/*             Surface Color Models             */
/************************************************/

/* This code is part of the surface library. */

#include "plsurf.h"
#include "plsurf_internal.h"

const plc_color  SURFACE_RED = {1,0,0,1}, 
  SURFACE_GREEN = {0,1,0,1}, 
    SURFACE_BLUE = {0,0,1,1},
      SURFACE_GRAY = {0.666,0.666,0.666,1}; 

void color_surface(surface *surf,
		   plc_color col)

     /* Procedure gives surface a constant color. */

{
  int i;

  for (i=0;i<surf->faces;i++) {

    surf->col_buf[i] = col;

  }

}

#ifdef CONVERTED

void color_by_proximity(surface *surf,
			polyline *curve,
			double radius)

     /* Procedure traces the polyline curve out on the surface */
     /* coloring faces near curve according to their "average" */
     /* distance from the curve. */

     /* The algorithm has two steps. In the first, we search for */
     /* edges near the curve, and check their lengths. If any is */
     /* longer than radius/5, it is subdivided, and the surface  */
     /* is retriangulated. */

     /* This process continues until all nearby edges are small  */
     /* enough to pass the test. */

     /* We then iterate through the faces of surf, coloring each by */
     /* proximity to the curve. Faces closest to the curve appear  */
     /* dark blue, while faces far from the curve are white. */

     /* We don't alter the color of any face at distance greater than */
     /* radius from the curve. */

{
  int *edge_list;
  int i,edges_added,j,good_dists;
  int small_enough = {FALSE};
  vector center;
  double distance,d[2];
  int divisions;

  /* Step 0: Alert the user. */

  fprintf(stdout,"color_by_proximity: Preparing to color surface."
                 " (%d verts,%d faces).\n",surf->verts,surf->faces);

  /* Step 1: Subdivide the nearby edges. */

  while (!small_enough) {
  
    if (surf->edges == NOT_COMPUTED) {

      make_edge_list(surf);

    }

    edge_list = calloc(surf->edges,sizeof(int));
    edges_added = 0;

    for(i=0;i<surf->edges;i++) {

      d[0] = distance_to_polyline(&surf->vert_buf[surf->edge_list[i][1]],
				  curve);

      d[1] = distance_to_polyline(&surf->vert_buf[surf->edge_list[i][2]],
				  curve);

      if (d[0] < radius && d[1] < radius) {

	/* do nothing. This is ok. */

      } else if (d[0] > radius && d[1] > radius) {

	/* do nothing. This is also ok. */
	
      } else if (edge_length(surf,i) > radius) {

	edge_list[edges_added++] = i;

      }
      
    }

    /* Now alert the user. */

    fprintf(stderr,"color_by_proximity: Preparing to subdivide %d edges.\n",
	    edges_added);

    /* Now go on. */

    if (edges_added > 0) {

      small_enough = FALSE;

      subdivide_edges(surf,edges_added,edge_list);      
      /* sanity_surface(surf); */

      fprintf(stderr,"color_by_proximity: Retriangulating surface.\n");

      triangulate_surface(surf);
      /* sanity_surface(surf); */

    } else {

      small_enough = TRUE;

    }

    free(edge_list);

  }
  
  /* Now alert the user again. */

  fprintf(stderr,"color_by_proximity: Dynamic subdivision complete. \n"
	         "                    New surface has %d verts, %d faces.\n",
	  surf->verts,surf->faces);

  /* Step 2: Now go through the list of faces, coloring by distance. */

  for(i=0;i<surf->faces;i++) {

    for(j=0,good_dists=0;j<3;j++) {

      if (distance_to_polyline(&surf->vert_buf[surf->face_buf[i][j]],curve) <
	  radius) {

	good_dists++;

      }

      if (good_dists == 3) {

	surf->col_buf[i].r = 0;
	surf->col_buf[i].g = 0;
	surf->col_buf[i].b = 0;
	surf->col_buf[i].a = 1.0;

      } else if (good_dists == 2) {

	surf->col_buf[i].r = 0.33;
	surf->col_buf[i].g = 0.33;
	surf->col_buf[i].b = 0.33;
	surf->col_buf[i].a = 1.0;

      } else if (good_dists == 1) {

	surf->col_buf[i].r = 0.66;
	surf->col_buf[i].g = 0.66;
	surf->col_buf[i].b = 0.66;
	surf->col_buf[i].a = 1.0;

      } else {

	surf->col_buf[i].r = 1;
	surf->col_buf[i].g = 1;
	surf->col_buf[i].b = 1;
	surf->col_buf[i].a = 1.0;

      }
  
    }

  }

}


#endif



