/**********************************************/
/*               Surface I/O                  */
/**********************************************/

/* This is part of the libsurf surface library. */

#include "plsurf.h"
#include "plsurf_internal.h"

int  load_surf_from_OFF(surface *surf,FILE *infile)

/*       Procedure reads a surface from a data file in the Geomview off */
/*       format. The procedure can deal with ordinary 3-dimensional off */
/*       files only. */

/*       There are a variety of color options available in an OFF file, */
/*       which makes reading them rather tricky. The procedure, as written, */
/*       will not handle entries into the Geomview colormap, but will  */
/*       correctly deal with integer, real, or missing color specifications. */

/*       The file infile is expected to be open, and ready for reading. */


{
  int test_char;
  int verts, faces, edges;
  int i,j,dim;

  plc_color SURFACE_BLUE = {0,0,1,1};
  plc_color SURFACE_GRAY = {0.6,0.6,0.6};
  plc_color thiscolor;

  /* First, we double-check that the file is open and legal. */

  if (infile == NULL) {

    fprintf(stderr,"load_surf_from_OFF: File not opened correctly.\n");
    return 0;

  }

  /* We check for the keyword "OFF". */

  test_char = fgetc(infile);

  if (test_char == 'O') {

    ungetc('O',infile);
    fscanf(infile,"OFF ");
  
  } else {

    ungetc(test_char,infile);

  }

  /* We could also have the keyword nOFF. This must be followed by the dimension. */

  test_char = fgetc(infile);

  if (test_char == 'n') {

    ungetc('n',infile);
    fscanf(infile,"nOFF %d",&dim);

    if (dim != 3) {

      fprintf(stderr,"load_surf_from_OFF: Cannot deal with nOFF of dimension %d != 3.\n",
	      dim);
      return 0;

    }

  } else {

    ungetc(test_char,infile);

  }

  /* Now, having taken care of the keyword, we */
  /* look for three integers, containing verts, faces, and edges. */

  if (scanints(infile,3,&verts,&faces,&edges) != 3) {

    fprintf(stderr,"load_surf_from_OFF: Cannot read line Verts Faces Edges"
	    " from file.\n");

    return 0;

  }

  if (verts <= 0 || faces <= 0) {

    fprintf(stderr,"load_surf_from_OFF: NVerts (%d) or NFaces (%d) bad.\n",
	    verts,faces);
    
    return 0;
  
  }

  /* We now initialize the surface.*/

  if (!new_surface(verts,faces,SURFACE_BLUE,surf)) {
  
    fprintf(stderr,"load_surf_from_OFF: Couldn't allocate surface.\n");
    return 0;

  }

  /* We now prepare to read NVerts triples of coordinates. */

  /* fprintf(stderr,"load_surf_from_OFF: Preparing to read %d vertices.\n",
     verts); */

  for(i=0;i<verts;i++) {

    if (scandoubles(infile,3,
	       &(surf->vert_buf[i].c[0]),
	       &(surf->vert_buf[i].c[1]),
	       &(surf->vert_buf[i].c[2])) != 3) {

      fprintf(stderr,"load_surf_from_OFF: Bad vertex (%d).\n",i);

      return 0;

    }

  }

  /*  fprintf(stderr,"load_surf_from_OFF: %d vertices read ok.\n",i); */

  /* Now we read the faces. */

  /*  fprintf(stderr,"load_surf_from_OFF: Preparing to read %d faces.\n",
      faces); */

  for(i=0;i<faces;i++) {

    /* First, we read the number of vertices on this face. */

    if (scanints(infile,1,&(surf->vtf[i])) != 1) {

      fprintf(stderr,"load_surf_from_OFF: Couldn't read number of vertices"
	      " in face %d.\n",i);

      return 0;

    }

    if (surf->vtf[i] <= 0) {

      fprintf(stderr,"load_surf_from_OFF: Face %d cannot contain %d"
	      " vertices.\n",i,surf->vtf[i]);

      return 0;

    }

    /* Now, we try to allocate a buffer that size. */

    surf->face_buf[i] = calloc(surf->vtf[i],sizeof(int));

    if (surf->face_buf[i] == NULL) {

      fprintf(stderr,"load_surf_from_OFF: Couldn't allocate buffer "
	      "containing %d vertices for face %d.\n",surf->vtf[i],i);

      return 0;

    }

    /* Having allocated the buffer, we are ready to read the numbers. */

    for(j=0;j<surf->vtf[i];j++) {

      if (scanints(infile,1,&(surf->face_buf[i][j])) != 1) {

	fprintf(stderr,"load_surf_from_OFF: Couldn't read vertex %d of "
		"face %d.\n",j,i);

	return 0;

      }

      if ((surf->face_buf[i][j] < 0) || (surf->face_buf[i][j] > (verts-1))) {

	fprintf(stderr,"load_surf_from_OFF: Vertex %d of face %d references"
		" illegal vertex number (%d).\n",j,i,surf->face_buf[i][j]);

	return 0;

      }

    }

    /* Now that the vertex references have been read, and passed, the */
    /* remaining challenge is the color information. For this, we are */
    /* sensitive to linebreaks in the OFF file. */
    
    if (!scancolor(infile,&thiscolor)) { /* Default color */
      
      surf->col_buf[i] = SURFACE_GRAY;
      
    } else {
      
      surf->col_buf[i] = thiscolor;
      
    }

    /* We are now done reading the face, and can close the loop. */

  }
      
  /* fprintf(stderr,"load_surf_from_OFF: %d faces read OK.\n",i); */

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

int color_compare(const void *a,const void *b)

{
  plc_color *A,*B;
  double coloreps = 0.01;

  A = (plc_color *)a;
  B = (plc_color *)b;

  if      (A->r < B->r - coloreps) { return -1; }
  else if (A->r > B->r + coloreps) { return 1; }
  else if (A->g < B->g - coloreps) { return -1; }
  else if (A->g > B->g + coloreps) { return 1; }
  else if (A->b < B->b - coloreps) { return -1; }
  else if (A->b > B->b + coloreps) { return 1; }
  else return 0;

}

plc_color *build_color_list(surface *surf,int *ncolors)

/* Builds a list of rgb colors used in a surface. Can be pretty slow. */

{
  plc_color *col_buf,*found;
  int nfound=0,i;
  
  col_buf = calloc(surf->faces,sizeof(plc_color));

  for(i=0;i<surf->faces;i++) {

    /* Check the list for this color */

    found = bsearch((void *)(&surf->col_buf[i]),
		    (void *)(col_buf),nfound,sizeof(plc_color),
		    color_compare);

    /* If not found, add it... */

    if (found == NULL) {

      assert(nfound < surf->faces);
      col_buf[nfound++] = surf->col_buf[i];
      qsort(col_buf,nfound,sizeof(plc_color),color_compare);

    }

  }

  /* We now have a list of colors */

  *ncolors = nfound;
  return col_buf;

}

void write_surf_to_POV(surface *surf,FILE *outfile,struct uvbuf *uvb)

     /* Procedure writes the surface <surf> to a POV-ray smooth triangle
	mesh. If uvb != NULL, include texture data in the mesh. */
{
  int i,j;
  plc_vector vnorm;

  /* First, we make sure that surf is triangulated, */
  /* and has edge and incidence lists. */

  if (outfile == NULL) {

    return;

  }

  if (!is_triangulated(surf)) {

    fprintf(stderr,"write_surf_to_POV: Only works for triangulated surfaces.");
    return;

  }

  /* Getting the color information to work is a bit of a pain. We need
     to assemble a list of colors from the faces, declare these as
     textures, and then refer to these declarations in the mesh code. */

  plc_color *col_buf;
  int ncolors;

  col_buf = build_color_list(surf,&ncolors);

  for(i=0;i<ncolors;i++) {

    fprintf(outfile,
	    "#declare OFFTex%d = texture { pigment { color rgb<%g,%g,%g> }}\n",
	    i,col_buf[i].r,col_buf[i].g,col_buf[i].b);

  }

  /* Now we are free to work. */

  fprintf(outfile,"\n"
	  "mesh {\n\n");

  for(i=0;i<surf->faces;i++) {

    fprintf(outfile,"  smooth_triangle {\n");

    for(j=0;j<3;j++) {

      fprintf(outfile,"    <%g,%g,%g>, ",
	      surf->vert_buf[surf->face_buf[i][j]].c[0],
	      surf->vert_buf[surf->face_buf[i][j]].c[1],
	      surf->vert_buf[surf->face_buf[i][j]].c[2] 
	      );
      
      vnorm = vert_normal(surf,surf->face_buf[i][j]);

      fprintf(outfile," <%g,%g,%g>",
	      vnorm.c[0],
	      vnorm.c[1],
	      vnorm.c[2]
	      );
      
      if (j != 2) {

	fprintf(outfile,",\n");

      }

    }

    if (uvb != NULL) {

      fprintf(outfile,"\n    uv_vectors <%g,%g>, <%g,%g>, <%g,%g>\n",
	      uvb->uv[surf->face_buf[i][0]].c[0],
	      uvb->uv[surf->face_buf[i][0]].c[1],
	      uvb->uv[surf->face_buf[i][1]].c[0],
	      uvb->uv[surf->face_buf[i][1]].c[1],
	      uvb->uv[surf->face_buf[i][2]].c[0],
	      uvb->uv[surf->face_buf[i][2]].c[1]);

    }
	      
    /* Now we add color information (if present in the plsurf object) */
    
    if (surf->col_buf != NULL) {

      int col_idx;
      plc_color *found;

      found = bsearch(&(surf->col_buf[i]),col_buf,ncolors,
		      sizeof(plc_color),color_compare);

      assert( found != NULL );
      col_idx = (int)(found - col_buf);

      fprintf(outfile,"\n    texture { OFFTex%d }\n",col_idx);

    }

    fprintf(outfile,"  }\n");

  }

  fprintf(outfile,"}\n");

}     


void print_edge_list(surface *surf)

  /* Procedure prints the edge list of the surface to stdout. */

{
  int i,j;

  if (surf->edges <= 0 && surf->edge_list == NULL) {

    printf("print_edge_list: Edge list not computed!\n");
    return;

  }

  printf("\nEdge List.\n\n");
  printf("# of faces  v_0  v_1  faces + orientations \n"
	 "----------  ---  ---  ---------------------------------\n");
  
  for(i=0;i<surf->edges;i++) {

    printf("%10d  %3d  %3d  ",
	   surf->edge_list[i][0],
	   surf->edge_list[i][1],
	   surf->edge_list[i][2]);

    for(j=0;j<surf->edge_list[i][0];j++) {

      printf("%3d  %2d  ",
	     surf->edge_list[i][2*j+3],
	     surf->edge_list[i][2*j+4]);

    }

    printf("\n");

  }

  printf("\n");

}

void print_incidence_lists(surface *surf)

  /* Procedure prints the incidence lists of the surface to stdout. */

{
  int i,j;

  if (surf->incident_faces == NULL || surf->incident_angles == NULL) {

    printf("print_incidence_lists: Lists not computed!\n");
    return;

  }

  printf("\nIncident Face and Angle Lists.\n\n");
  printf("# of faces  faces or angles \n"
	 "----------  ---------------------------------\n");
  
  for(i=0;i<surf->verts;i++) {

    printf("%10d ",
	   surf->ftv[i]);

    for(j=0;j<surf->ftv[i];j++) {

      if (surf->incident_faces[i][j] != NO_FACE) {
	
	printf("%6d  ",
	       surf->incident_faces[i][j]);

      } else {

	printf(" NO_FACE ");

      }

    }

    printf("\n");

    printf("%10d ",
	   surf->ftv[i]);
    
    for(j=0;j<surf->ftv[i];j++) {

      printf("%6.3f  ",
	     surf->incident_angles[i][j]);
      
    }

    printf("\n\n");

  }
  
  printf("\n");

}

#ifdef CONVERTED

void print_shared_faces(surface *surf)

     /* Procedure prints the shared faces of every pair of edges */
     /* vertices, and vertex/edges of the surface to stdout. */

{
  int i,j,k;
  int num,*result;

  if (surf->incident_faces == NULL || surf->incident_angles == NULL ||
      surf->edge_list == NULL) {

    printf("print_shared_faces: Lists not computed!\n");
    return;

  }

  printf("\nShared Face Lists.\n\n");

  printf("(edge,edge)  faces \n"
	 "-----------  ---------------------------------\n");
  
  for(i=1;i<surf->edges;i++) {

    for(j=0;j<i;j++) {

      printf("(%4d,%4d) ",
	     i,j);

      num = edge_shared_faces(surf,i,j,&result);

      for(k=0;k<num;k++) {

	printf("%4d  ",result[k]);

      }

      free(result);

      printf("\n");

    }

  }
  
  printf("\n");

  printf("(vert,vert)  faces \n"
	 "-----------  ---------------------------------\n");
  
  for(i=1;i<surf->verts;i++) {

    for(j=0;j<i;j++) {

      printf("(%4d,%4d) ",
	     i,j);

      num = vert_shared_faces(surf,i,j,&result);

      for(k=0;k<num;k++) {

	printf("%4d  ",result[k]);

      }

      free(result);

      printf("\n");

    }

  }
  
  printf("\n");

  printf("(edge,vert)  faces \n"
	 "-----------  ---------------------------------\n");
  
  for(i=1;i<surf->edges;i++) {

    for(j=0;j<surf->verts;j++) {

      printf("(%4d,%4d) ",
	     i,j);

      num = edge_vert_shared_faces(surf,i,j,&result);

      for(k=0;k<num;k++) {

	printf("%4d  ",result[k]);

      }

      free(result);

      printf("\n");

    }

  }
  
  printf("\n");

}

#endif
