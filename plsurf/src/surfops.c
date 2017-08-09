/**********************************************/
/*            Surface Operations              */
/**********************************************/

/* This code is part of the surface library. */

#include"plsurf.h"
#include"plsurf_internal.h"

void scale_surface(surface *surf,double scale) 

     /* Procedure scales surf by <scale>. */

{
  int i;

  for(i=0;i<surf->verts;i++) {

    surf->vert_buf[i] = plc_scale_vect(scale,surf->vert_buf[i]);

  }

}

#ifdef CONVERTED 

void normalize_surface_volume(surface *surf)

     /* Procedure scales surface to unit volume. */

{ 
  scale_surface(surf,1/pow(volume_surface(surf),1.0/3.0));
}

#endif
		
void translate_surface(surface *surf,plc_vector vec)
     
     /* Procedure moves a surface by vec. */

{
  int i;

  for(i=0;i<surf->verts;i++) {

    surf->vert_buf[i] = plc_vect_sum(surf->vert_buf[i],vec);

  }
}


void spherical_projection_surface(surface *surf)

     /* This procedure moves all the vertices of surf to their 
	projected positions on the unit sphere centered at the
	origin. */

{
  int  i;
  bool ok;
 
  for(i=0;i<surf->verts;i++) {

    surf->vert_buf[i] = plc_normalize_vect(surf->vert_buf[i],&ok);

    if (!ok) {

      fprintf(stderr,"spherical_projection_surface: Vertex %d too "
	      "close to origin to project.",i);
      exit(1);

    }
  }

}

#ifdef CONVERTED

void   unroll_cone(surface *cone)

     /* We assume that cone has been produced by make_cone, from the polyline
	library. Using some information about how this procedure works, we can
	compute the development of this surface onto the plane. 

	This only works if the cone is open. To produce an open cone from a closed
	polyline, set the pline->closed flag to FALSE before coning. */

{
  int n_rays;
  int i,ray,sf_cnt;
  
  int    *ray_verts;
  double *ray_angles,ra_so_far,radius;

  surfpt vertex;
  int    closed;

  /* First, we make sure the cone is open */

  make_edge_list(cone);
  make_incidence_lists(cone);
  vertex.type = VERTPT; vertex.vert = 0; vertex.surf = cone;

  if (!on_bdy(vertex)) {

    fprintf(stderr,"unroll_cone: Can't unroll a closed cone.\n");
    exit(1);

  }
  
  /* Next, we count the total number of rays, keeping in mind that each face
     involving vertex 0 indicates a ray. */

  for(i=0,n_rays=2;i<cone->faces;i++) {

    if (cone->face_buf[i][0] == 0) {n_rays++;}

  }

  /* We have now counted rays. We can now figure out the angles and vertices. */

  ray_verts = calloc(n_rays,sizeof(int));
  ray_angles = calloc(n_rays,sizeof(double));
  ra_so_far  = 0.0;

  ray_verts[0]  = 1;		/* We set up the first ray. It starts immediately at angle 0. */
  ray_angles[0] = 0.0;
  sf_cnt        = 1;

  ray_verts[n_rays-1]  = cone->verts+1; /* We make sure that we never reach this. */
  ray_angles[n_rays-1] = 0;

  for(i=0;i<cone->faces;i++) {

    if (cone->face_buf[i][0] == 0) { /* We are at a special face. */

      ray_verts[sf_cnt]  = cone->face_buf[i][2];

      ra_so_far += face_angle_at_vert(cone,i,0);
      ray_angles[sf_cnt] = ra_so_far;

      sf_cnt++;

    }

  }

  /* We have now prepared a table of index vertices and angles. */
  /* It is easy to assign each vertex it's new location. */

  for(i=1,ray=0;i<cone->verts;i++) {

    if (i == ray_verts[ray+1]) {	/* Check to see if we've changed rays. */

      ray++;

    }

    radius = segmentlength(&cone->vert_buf[i],&cone->vert_buf[0]);
    
    cone->vert_buf[i].x = radius * cos(ray_angles[ray]);
    cone->vert_buf[i].y = radius * sin(ray_angles[ray]);
    cone->vert_buf[i].z = 0.0;

  }

  /* Now we can fix the 0th vertex. */

  cone->vert_buf[0].x = cone->vert_buf[0].y = cone->vert_buf[0].z = 0.0;

  /* We are now done. */

  free(ray_verts); free(ray_angles);

}

surface join_surfaces(int n_surfs, ... )

     /* Procedure joins a collection of surfaces, returning the result in a new surface. */
     /* We use the va_arg syntax, so the list of pointers to surfaces should start after */
     /* the number n_surfs of surfaces to join. */
{
  int i,j,k;
  va_list ap;
  surface **surfs,ret;
  int verts = {0},faces = {0},v_ofs = {0},f_ofs={0};
  
  /* First, we do some error checking. */

  if (n_surfs < 1) {

    fprintf(stderr,"join_surfaces: Can't join %d surfaces.\n",n_surfs);
    exit(2);

  }

  va_start(ap,n_surfs);

  /* Now we load up the arguments. */

  surfs = calloc(n_surfs,sizeof(surface *));
  for(i=0;i<n_surfs;i++) {

    surfs[i] = va_arg(ap,surface *);

  }
  va_end(ap);

  /* We are prepared to work. */
  
  for(i=0;i<n_surfs;i++) {

    verts += surfs[i]->verts;
    faces += surfs[i]->faces;

  }

  new_surface(verts,faces,default_color(),&ret);

  for(i=0;i<n_surfs;i++) {

    for(j=0;j<surfs[i]->verts;j++) {

      ret.vert_buf[j+v_ofs] = surfs[i]->vert_buf[j];
      
    }

    for(j=0;j<surfs[i]->faces;j++) {

      ret.vtf[j+f_ofs] = surfs[i]->vtf[j];
      ret.face_buf[j+f_ofs] = calloc(ret.vtf[j+f_ofs],sizeof(int));

      for(k=0;k<ret.vtf[j+f_ofs];k++) {

	ret.face_buf[j+f_ofs][k] = surfs[i]->face_buf[j][k] + v_ofs;

      }

    }

    f_ofs += surfs[i]->faces;
    v_ofs += surfs[i]->verts;

  }

  free(surfs);
  return ret;

}
      
surface pushoff_surface(surface *surf,double eps)

     /* Procedure pushes each vertex in the direction of its vertex normal, */
     /* and makes a copy of the resulting surface. Computed structures such */
     /* as edge and incidence lists are not copied. */

{
  int i,j;
  surface ret;
  vector  vnorm;

  new_surface(surf->verts,surf->faces,default_color(),&ret);

  for(i=0;i<surf->verts;i++) {

    vnorm = vert_normal(surf,i);
    linear_combine(1.0,&surf->vert_buf[i],eps,&vnorm,&ret.vert_buf[i]);
    
  }

  for(i=0;i<surf->faces;i++) {

    ret.vtf[i] = surf->vtf[i];
    ret.face_buf[i] = calloc(ret.vtf[i],sizeof(int));

    for(j=0;j<ret.vtf[i];j++) {

      ret.face_buf[i][j] = surf->face_buf[i][j];

    }

  }

  return ret;

}
  
void    reverse_surface(surface *surf)

     /* Procedure reverses the orientation of a surface. */

{
  int i;
  int j;
  int *new_buf;

  for(i=0;i<surf->faces;i++) {

    new_buf = calloc(surf->vtf[i],sizeof(int));

    for(j=0;j<surf->vtf[i];j++) {

      new_buf[surf->vtf[i]-1-j] = surf->face_buf[i][j];

    }

    free(surf->face_buf[i]);
    surf->face_buf[i] = new_buf;

  }

}

surface *unfold_surface(surface *surf,int face,int edge)

     /* Procedure attempts an "orange peel" unfolding of a polyhedral surface
	starting with the given face and edge. The result is a planar surface. */

{
  int *visited,i;

  /* First, we check the stupid stuff. */

  if (surf->edges == NOT_COMPUTED) {

    make_edge_list(surf);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"unfold_surface: Face <%d> is out of range.\n",face);
    exit(1);

  }

  if (!edge_legal(surf,edge)) {

    fprintf(stderr,"unfold_surface: Edge <%d> is out of range.\n",edge);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"unfold_surface: Edge <%d> is not on face <%d>.\n",edge,face);
    exit(1);

  }

  /* Now we are ready to get to work. We begin with some initializations. */

  visited = calloc(surf->faces,sizeof(int));

  if (visited == NULL) {

    fprintf(stderr,"unfold_surface: Couldn't allocate 'visited' buffer.");
    exit(1);

  }

}  
  
void     close_surface_by_cone(surface *surf,vector vertex)

     /* Procedure closes the surface by coning each boundary component to vertex. */
     /* The cones constructed are quite simple. We attempt to respect the orientation */
     /* on the surface, if one exists. The procedure is destructive to surf. */

{
  vector *new_vert_buf;

  int  **new_face_buf;
  int   *new_vtf;
  int    i,j;
  color *new_col_buf;

  int  *bdy_edges,nbe;

  nbe = get_boundary_edges(surf,&bdy_edges);
  
  if (nbe == 0) {

    return;

  }

  /* Now we allocate new, bigger buffers for all the new faces. */

  new_vtf      = calloc(nbe+surf->faces,sizeof(int));
  new_col_buf  = calloc(nbe+surf->faces,sizeof(color));
  new_face_buf = calloc(nbe+surf->faces,sizeof(int *));

  for(i=0;i<surf->faces;i++) {

    new_vtf[i] = surf->vtf[i];
    new_col_buf[i] = surf->col_buf[i];
    
    new_face_buf[i] = calloc(new_vtf[i],sizeof(int));
    
    for(j=0;j<surf->vtf[i];j++) {

      new_face_buf[i][j] = surf->face_buf[i][j];

    }

  }

  new_vert_buf = calloc(surf->verts+1,sizeof(vector));

  for(i=0;i<surf->verts;i++) {

    new_vert_buf[i] = surf->vert_buf[i];

  }

  /* We now fill in the new vertex. */

  new_vert_buf[surf->verts] = vertex;

  /* We now fill in the new faces. */

  for(i=0;i<nbe;i++) {

    new_vtf[i+surf->faces] = 3;
    new_face_buf[i+surf->faces] = calloc(3,sizeof(int));

    if (surf->edge_list[bdy_edges[i]][4] == 1) { /* Positively oriented on the old face */
  
      new_face_buf[i+surf->faces][0] = surf->edge_list[bdy_edges[i]][2];
      new_face_buf[i+surf->faces][1] = surf->edge_list[bdy_edges[i]][1];

    } else {

       new_face_buf[i+surf->faces][0] = surf->edge_list[bdy_edges[i]][1];
       new_face_buf[i+surf->faces][1] = surf->edge_list[bdy_edges[i]][2];

    }

    new_face_buf[i+surf->faces][2] = surf->verts; /* The third is always the new vertex. */

    /* We now fix the color buffer to be the same as the color of the adjacent face. */

    new_col_buf[i+surf->faces] = surf->col_buf[surf->edge_list[bdy_edges[i]][3]];
    
  }

  /* We are now done. All that remains is to swap in the new buffers and free the old ones. */

  free(surf->vert_buf);
  free(surf->col_buf);
  free(surf->vtf);
  free(surf->face_buf);

  surf->verts = surf->verts + 1;
  surf->faces = surf->faces + nbe;

  surf->vert_buf = new_vert_buf;
  surf->vtf = new_vtf;
  surf->face_buf = new_face_buf;
  surf->col_buf = new_col_buf;

  surf->is_closed = TRUE;
  
}
  
#endif    
