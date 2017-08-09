/**********************************************/
/*            Standard Surfaces               */
/**********************************************/

/* This code is part of the surface library. */

#include"plsurf.h"
#include"plsurf_internal.h"


surface icosahedron()

     /* Procedure returns a surface containing a standard icosahedron. */

{
  surface surf;
  int     i,j;

  plc_vector  vert_buf[12] = { 
    {{0.0, 0.0, 2.0}},
    {{1.788854, 0.000000, 0.894427}},
    {{0.552786, 1.701302, 0.894427}},
    {{-1.447214, 1.051462, 0.894427}},
    {{-1.447214, -1.051462, 0.894427}},
    {{0.552786, -1.701302, 0.894427}},
    {{1.447214, 1.051462, -0.894427}},
    {{-0.552786, 1.701302, -0.894427}},
    {{-1.788854, 0.000000, -0.894427}},
    {{-0.552786, -1.701302, -0.894427}},
    {{1.447214, -1.051462, -0.894427}},
    {{0.0, 0.0, -2.0}}};

  int    vtf[20] = {3,3,3,3,3,
		    3,3,3,3,3,
		    3,3,3,3,3,
		    3,3,3,3,3};

  int    face_buf[20][3] = {{2,0,1},{3,0,2},{4,0,3},{5,0,4},{1,0,5},
			    {2,1,6},{7,2,6},{3,2,7},{8,3,7},{4,3,8},
			    {9,4,8},{5,4,9},{10,5,9},{6,1,10},{1,5,10},
			    {6,11,7},{7,11,8},{8,11,9},{9,11,10},{10,11,6}};

  new_surface(12,20,SURFACE_BLUE,&surf);

  for(i=0;i<surf.verts;i++) {
    
    surf.vert_buf[i] = vert_buf[i];

  }

  for(i=0;i<surf.faces;i++) {

    surf.vtf[i] = vtf[i];
    surf.face_buf[i] = calloc(surf.vtf[i],sizeof(int));
    
    for(j=0;j<surf.vtf[i];j++) {

      surf.face_buf[i][j] = face_buf[i][j];

    }

  }

  return surf;
}
	      
surface geodesic_dome(int subdivision,int recursion)

     /* Procedure creates a geodesic dome with 20*4*subdivision*recursion vertices
	by subdividing the faces of an icosahedron <subdivision> times, projecting
	them to the sphere, and repeating <recursion> times. */

{
  int i,j;

  surface dome;

  /* We first check the inputs. */

  if (subdivision < 0 || recursion < 0) {

    fprintf(stderr,"geodesic_dome: Illegal subdivision %d or recursion %d.\n",
	    subdivision,recursion);

    exit(2);

  }

  dome = icosahedron();

  for(i=0;i<recursion;i++) {

    for(j=0;j<subdivision;j++) {

      subdivide_faces(&dome);

    }

    spherical_projection_surface(&dome);

  }

  return dome;

}
