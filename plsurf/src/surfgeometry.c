/***********************************************/
/*       Computing Surface Geometry            */
/***********************************************/

/* This code is part of the surfaces library.  */

#include "plsurf.h"
#include "plsurf_internal.h"

#ifdef CONVERTED 

int inside(surface *surf,vector *point)

  /* Procedure decides whether point is inside the surface surf. */
  
  /* The algorithm works by computing the projection of surf to  */
  /* the sphere surrounding point using gauss_surface. If the    */
  /* integral is equal to 1, the point is inside, if the integral*/
  /* is equal to 0, the point is outside. */

  /* If the surface has any component with boundary, the results    */
  /* of this computation are unreliable. We check is_closed, before */
  /* making the computation. */

  /* We have written the code so that it does not depend on the  */
  /* orientation of the surface, though the surface must be oriented */
  /* in order to make the computation. */

{
  double result;

  if (!is_closed(surf)) {

    fprintf(stderr,"inside: Surface is not closed!\n");
    return FALSE;

  }

  if (!is_oriented(surf)) {

    fprintf(stderr,"inside: Surface is not oriented!\n");
    return FALSE;

  }

  result = gauss_surface(surf,point);

  if (fabs(result) < 1.1 && fabs(result) > 0.9) {

    return TRUE;
    
  } else if (fabs(result) < 0.1) {

    return FALSE;

  } else {

    fprintf(stderr,"inside: Forbidden gauss integral value of %g. \n",
	    result);

    return FALSE;

  }

}

double gauss_surface(surface *surf,vector *origin)

  /* Procedure computes the area of the projection of surf onto the */
  /* sphere centered at origin, under the Gauss map. This depends on */
  /* the orientation of the surface, so we require that the surface */
  /* carry a consistent orientation to make the computation. */

  /* We divide the answer by 4 Pi, to correspond to the area of the unit */
  /* sphere. */

{
  vector *temp_verts;
  int    face_cnt,vert_cnt;
  double running_sum = {0};

  if (!is_oriented(surf)) {

    fprintf(stderr,"gauss_surface: Surface is not correctly oriented! \n"
	           "               Cannot compute Gauss integral.\n");

    fprintf(stderr,"gauss_surface: Terminating program.\n");

    exit(1);

  }

  for(face_cnt=0;face_cnt<surf->faces;face_cnt++) {

    temp_verts = calloc(surf->vtf[face_cnt],sizeof(vector));

    for(vert_cnt=0;vert_cnt<surf->vtf[face_cnt];vert_cnt++) {

      temp_verts[vert_cnt] = 
	surf->vert_buf[surf->face_buf[face_cnt][vert_cnt]];
    
    }

    running_sum += polygon_proj(origin,surf->vtf[face_cnt],temp_verts);

    free(temp_verts);

  }

  return (running_sum/(2*TWO_PI));

}

double polygon_proj(vector *org,     /* origin to project on */ 
		    int n,           /* # of sides of polygon */
		    vector *verts    /* array of vertices */
		    ) 

  /* Computes the area of the projection of a planar polygon onto a unit 
     sphere centered at org. The function assumes that the vertices are
     given in counterclockwise order. If the polygon is "above" org, the
     output is positive, and if it's "below", the output is negative. */

{

  double area = 0;         /* running sum for area of projection */

  int i;

  for (i=1; i<n-1 ;i++) {

    area += triangle_proj(org, &verts[0], &verts[i], &verts[i+1]);
  
  }

  return area;
}


double triangle_proj(vector *org,     /* origin to project on */ 
		     vector *a,       /* first point */
		     vector *b,       /* second point */
		     vector *c        /* third point */
		     ) 

  /* Computes the area of the projection of a triangle onto a unit sphere
     centered at org. The function assumes that the vertices are given in
     counterclockwise order. If the triangle is "above" org, the output is 
     positive, and if it's "below", the output is negative. */

{
  vector normal;          /* normal vector to the triangle, determined by
			     the orientation of the boundary */
  vector pts[4];          /* the points */
  vector diffs[4][4];     /* differences between the points */
  vector n[2];            /* normal vectors to planes */
  double area = {-PI};    /* running sum for area of projection */
  double predicted_area;  /* area prediction by altitudes to org */
  double lengths[3];      /* side lengths of the triangle; lengths[i] is the
			     length of the side from pts[i-1] to pts[i] */
  double alts[3];         /* altitudes from the sides to the origin */
  double sign;            /* whether the result should be + or - */
  int i,j;
  double DELTA = {0.0000000001};
  
  pts[0] = *a, pts[1] = *b, pts[2] = *c, pts[3] = *org;
  
  for(i=0;i<4;i++) {
    
    for(j=0;j<4;j++) {
    
      vectordiff(&pts[i],&pts[j],&diffs[i][j]);
    
    }
    
  }

  for(i=0;i<3;i++) {
    
    j = i + 3;
    lengths[i] = norm(&diffs[j%3][(j-1)%3]);

  }

  cross(&diffs[1][0], &diffs[2][1], &normal);
  
  /* error trapping */

  if (fabs(dot(&diffs[0][3], &normal)) < DELTA) {


    /* if org is almost in the same plane as the triangle */
    
    /* check if it's the same as one of the endpoints */
    for(i=0; i<3; i++) {
      if(segmentlength(&pts[i], &pts[3]) < DELTA)
	return 0;
    }
    
    /* otherwise, see if it's inside the triangle */
    area = triangle_area(a, b, c);
    for(i=0; i<3; i++) {
      
      j = i + 3;
      alts[i] = dist_to_line(org, &pts[(j-1)%3], &pts[j%3]);
      
    }

    predicted_area = (lengths[0]*alts[0] + lengths[1]*alts[1] + 
		      lengths[2]*alts[2]) / 2;
    
    if(fabs(predicted_area - area) < DELTA) {
      
      /* org is inside the triangle; print out an error message */

      fprintf(stderr, "triangle_proj: Could not calculate the projection.\n"
	              "               Origin is inside triangle.\n");
      
    }
    
    /* If org is outside the triangle, the projection is truly 0. Otherwise,
       we return 0 after printing the error message. */
    
    return 0;
  }
  
  /* figure out the sign of the result*/
  if(dot(&diffs[0][3], &normal) > 0)
    sign = 1;
  else
    sign = -1;

  /* now, we can finally calculate the area */
  
  for (i=0;i<3;i++) {    

    j = i + 3;
    cross(&diffs[(j-1)%3][3], &diffs[j%3][3], &n[0]);
    scalarmultiply(sign, &n[0]);

    /* this ensures that we are working with inward-pointing normals */

    cross(&diffs[j%3][3], &diffs[(j+1)%3][3], &n[1]);
    scalarmultiply(sign, &n[1]);
    area += PI - acos(dot(&n[0],&n[1])/(norm(&n[0])*norm(&n[1])));
  }
  
  return(area * sign);
}      


double area_surface(surface *surf)

  /* Procedure computes the area of a surface. */

{
  int i;
  double run_sum = {0};

  for (i=0;i<surf->faces;i++) {

    run_sum += area_face(surf,i);

  }

  surf->area = run_sum;

  return run_sum;

}

double volume_surface(surface *surf) 

     /* Computes the volume enclosed by a closed surface. */
     /* Requires that the surface be triangulated. */

{
  double vol = {0};
  int    i;

  if (!is_triangulated(surf)) {

    fprintf(stderr,"volume_surface: Surface must be triangulated.\n");
    exit(1);

  }

  if (!is_oriented(surf)) {

    fprintf(stderr,"volume_surface: Surface must be oriented.\n");
    exit(1);

  }

  if (!is_closed(surf)) {

    fprintf(stderr,"volume_surface: Surface must be closed.\n");
    exit(1);

  }

  /* Now we are ready to work. */

  for(i=0;i<surf->faces;i++) {

    vol += triple(&surf->vert_buf[surf->face_buf[i][0]],
		  &surf->vert_buf[surf->face_buf[i][1]],
		  &surf->vert_buf[surf->face_buf[i][2]]) / 2.0;

  }

  return fabs(vol);

}

double area_face(surface *surf,int i)

  /* This procedure computes the area of a (planar) face of a surface. */
  /* The algorithm works by generating a new frame, in which the first */
  /* edge of the polygon lies on the positive x-axis, and the surface  */
  /* normal is the z-axis. */

  /* We then change basis for all the vectors in the face, and compute */
  /* the area using the computational geometry trick of computing the  */
  /* flux of the vector field x (x_hat) (which has unit divergence)    */
  /* across the perimeter. */

{
  matrix basis;
  vector *diffs,*new_vects;
  int n,j;
  double total = {0};
  vector triangle[3];

  if (surf->vtf[i] == 3) {

    triangle[0] = surf->vert_buf[surf->face_buf[i][0]];
    triangle[1] = surf->vert_buf[surf->face_buf[i][1]];
    triangle[2] = surf->vert_buf[surf->face_buf[i][2]];

    return triangle_area(&triangle[0],&triangle[1],&triangle[2]);

  }

  if (!face_planar(surf,i)) {

    fprintf(stderr,"area_face: Cannot compute area for face %d\n",
	           "           It does not lie in a single plane\n",i);
    exit(1);

  }

  /* Our first step is to move all the vectors in question so that the  */
  /* first vertex lies at the origin. */

  diffs = calloc(surf->vtf[i],sizeof(vector));
  new_vects = calloc(surf->vtf[i],sizeof(vector));

  if (diffs == NULL || new_vects == NULL) {

    fprintf(stderr,"area_face: Could not allocate %d vertex buffer.\n",
	           "           Terminating.\n",
	    surf->vtf[i]);

    exit(1);

  }

  for(j=0;j<surf->vtf[i];j++) {

    vectordiff(&surf->vert_buf[surf->face_buf[i][j]],
	       &surf->vert_buf[surf->face_buf[i][0]],
	       &diffs[j]);

  }

  /* Now, we must generate the basis matrix. */

  basis.rows[0] = diffs[1];
  basis.rows[2] = face_normal(surf,i);
  cross(&basis.rows[2],&basis.rows[0],&basis.rows[1]);

  for(j=0;j<3;j++) {
   
    vectornormalize(&basis.rows[j]);
  
  }

  /* We can now rewrite each vector in the new basis. */

  for(j=0;j<surf->vtf[i];j++) {

    change_basis(&diffs[j],&basis,&new_vects[j]);

  }

  /* Last, we'll compute the area. */

  for(j=0;j<surf->vtf[i]-1;j++) {

    total += (0.5)*(new_vects[j].x * new_vects[j+1].y) 
      + (0.5)*(new_vects[j+1].x * new_vects[j+1].y) 
      - (0.5)*(new_vects[j].x * new_vects[j].y)
      - (0.5)*(new_vects[j+1].x * new_vects[j].y);
    
  }

  total +=  (0.5)*(new_vects[j].x * new_vects[0].y) 
      + (0.5)*(new_vects[0].x * new_vects[0].y) 
      - (0.5)*(new_vects[j].x * new_vects[j].y)
      - (0.5)*(new_vects[0].x * new_vects[j].y);

  return total;

}

double edge_length(surface *surf,int edge) 

     /* Procedure returns the length of the ith edge of surface. */

{
  /* First, we check the edge_list. */

  if (surf->edge_list == NULL) {

    fprintf(stderr,"edge_length: Edge list not allocated.\n");
    return -1.0;

  }

  /* Now we compute. */

  return segmentlength(&(surf->vert_buf[surf->edge_list[edge][1]]),
		       &(surf->vert_buf[surf->edge_list[edge][2]]));

}

void   face_outline(surface *surf,int i,polyline *pline) 

     /* Procedure writes the outline of the face to a polyline. */

{
  int j;
  vector OUTLINE = {1,1,1};

  new_polyline(pline,surf->vtf[i]+1,"outline",OUTLINE);

  for(j=0;j<surf->vtf[i];j++) {

    pline->buffer[j] = surf->vert_buf[surf->face_buf[i][j]];

  }

  pline->buffer[j] = surf->vert_buf[surf->face_buf[i][0]];
  closed_polyline(pline);

}

#endif

plc_vector face_normal(surface *surf,int i)

  /* Procedure computes a normal vector to the ith face of surf, */
  /* assuming that this face is planar. The code is smart enough */
  /* to handle faces where edges lie in a single line. */

{
  int j;
  plc_vector normal,v[3];
  int vts;

  /* First, we check the stupid stuff. */

  if (!face_legal(surf,i)) {

    fprintf(stderr,"face_normal: Illegal face index %d.\n",i);
    exit(1);

  }

  /* Now we work. */

  vts = surf->vtf[i];

  for(j=0;j<surf->vtf[i];j++) {

    v[0] = surf->vert_buf[surf->face_buf[i][j]];
    v[1] = surf->vert_buf[surf->face_buf[i][(j+1)%vts]];
    v[2] = surf->vert_buf[surf->face_buf[i][(j+2)%vts]];
    
    bool ok;
    
    normal = plc_normal(v[0],v[1],v[2],&ok);

    if (ok) { return normal; }
 
  }

  fprintf(stderr,"face_normal: Error! Can't compute normal vector on face %d.\n",i);
  exit(1);

} 

#ifdef CONVERTED 

vector edge_normal(surface *surf,int i)

     /* Procedure computes an interpolated normal to the surface */
     /* along an edge. Of course, we need to have computed the */
     /* edge_list for the index i to make sense. */

{
  vector sum = {0,0,0},this_normal;
  int j;

  /* First, we check to make sure that the edge list has been made. */
  
  if (surf->edge_list == NULL) {

    fprintf(stderr,"edge_normal: Edge list not computed!\n");
    return;

  }

  /* Now we find the face normals of the faces on this edge, and add them */

  for(j=0;j<surf->edge_list[i][0];j++) {

    this_normal = face_normal(surf,surf->edge_list[i][2*j + 3]);
    vectoradd(&sum,&this_normal);

  }

  vectornormalize(&sum);

  return sum;

}

#endif

plc_vector vert_normal(surface *surf,int i)

     /* Procedure requires that we've computed the incidence lists. */
     /* If so, it returns an interpolated normal vector at this vertex. */

{
  int j;
  plc_vector this_normal,sum={{0,0,0}};

  /* First, we check that the incidence lists have been made. */

  if (surf->incident_faces == NULL) {

    make_incidence_lists(surf);

  }

  /* Next, we gather the normals from the faces incident to i, and sum. */

  for(j=0;j<surf->ftv[i];j++) {

    if (surf->incident_faces[i][j] != NO_FACE) {

      this_normal = face_normal(surf,surf->incident_faces[i][j]);
      sum = plc_vect_sum(sum,this_normal);

    }

  }

  if (surf->ftv[i] == 0) {	/* This could be a "stray" vertex. */
                                /* In this case, assign it (0,0,1). */

    sum = plc_build_vect(0,0,1);

  } else {

    bool ok;
    sum = plc_normalize_vect(sum,&ok);
    if (!ok) {fprintf(stderr,"vert_normal: Can't compute vertex normal at vert %d.",i); exit(1); }

  }

  return sum;

}

#ifdef CONVERTED

vector inward_normal(surface *surf,int edge,int face)

     /* Procedure returns the (unit) inward pointing normal of
	EDGE as an edge of FACE, in the boundary orientation of
	EDGE on FACE.

	If EDGE does not bound FACE, we return an error and quit. */

{
  vector edge_vec,f_normal,result;
  double sign;
  
  /* As always, first we check the dumb stuff. */

  if (surf == NULL) {

    fprintf(stderr,"inward_normal: NULL surface passed.\n");
    exit(1);

  }

  if (!edge_legal(surf,edge)) {

    fprintf(stderr,"inward_normal: Edge number %d not legal.\n",edge);
    exit(1);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"inward_normal: Face number %d not legal.\n",face);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"inward_normal: Edge %d not on face %d.\n",edge,face);
    exit(1);

  }

  /* Now we are ready to work. We create a vector pointing along the edge. */

  sign = (double)(edge_orientation_on_face(surf,edge,face));

  linear_combine(-1*sign,&surf->vert_buf[surf->edge_list[edge][1]],
		  1*sign,&surf->vert_buf[surf->edge_list[edge][2]],
		 &edge_vec);

  vectornormalize(&edge_vec);

  f_normal = face_normal(surf,face);

  /* The inward normal is now the cross product of f_normal and edge_vec. */
  
  cross(&f_normal,&edge_vec,&result);

  /* We are now done. */

  return result;

}

vector outward_normal(surface *surf,int edge,int face)

     /* Procedure returns the (unit) outward pointing normal of
	EDGE as an edge of FACE, in the boundary orientation of
	EDGE on FACE.

	If EDGE does not bound FACE, we return an error and quit. */

{
  vector edge_vec,f_normal,result;
  double sign;
  
  /* As always, first we check the dumb stuff. */

  if (surf == NULL) {

    fprintf(stderr,"outward_normal: NULL surface passed.\n");
    exit(1);

  }

  if (!edge_legal(surf,edge)) {

    fprintf(stderr,"outward_normal: Edge number %d not legal.\n",edge);
    exit(1);

  }

  if (!face_legal(surf,face)) {

    fprintf(stderr,"outward_normal: Face number %d not legal.\n",face);
    exit(1);

  }

  if (!edge_on_face(surf,edge,face)) {

    fprintf(stderr,"outward_normal: Edge %d not on face %d.\n",edge,face);
    exit(1);

  }

  /* Now we are ready to work. We create a vector pointing along the edge. */

  sign = (double)(edge_orientation_on_face(surf,edge,face));

  linear_combine(-1*sign,&surf->vert_buf[surf->edge_list[edge][1]],
		  1*sign,&surf->vert_buf[surf->edge_list[edge][2]],
		 &edge_vec);

  vectornormalize(&edge_vec);

  f_normal = face_normal(surf,face);

  /* The inward normal is now the cross product of edge_vec and f_normal. */
  
  cross(&edge_vec,&f_normal,&result);

  /* We are now done. */

  return result;

}

vector barycenter(surface *surf,int face)

     /* Returns the barycenter to this face: that is, the average position
	of all the vertices. */

{
  double weight;
  vector run_sum = {0,0,0};
  int i;

  weight = 1.0/(double)(surf->vtf[face]);

  for(i=0;i<surf->vtf[face];i++) {

    linear_combine(1.0,&run_sum,
		   weight,&surf->vert_buf[surf->face_buf[face][i]],
		   &run_sum);

  }

  return run_sum;
}

#endif

plc_vector vertex_center_of_mass(surface *surf)

     /* Procedure computes the center of mass of the vertices of surf. */

{
  int i;
  plc_vector com = {{0,0,0}};

  if (surf == NULL) {

    fprintf(stderr,"vertex_center_of_mass: Passed bad surface data.\n");
    exit(2);

  }

  for(i=0;i<surf->verts;i++) {

    com = plc_vect_sum(com,surf->vert_buf[i]);
    
  }

  com = plc_scale_vect(1.0/(double)surf->verts,com);
  return com;
}
  

void bbox_surface(surface *surf,plc_vector *bottom,plc_vector *top)

  /* Procedure determines a bounding box for the surface. */
  /* Since all the faces are planar, this merely involves */
  /* searching the space of vertices for their least and  */
  /* greatest coordinates. */

     /* We must be wary of the case where some dimension of this */
     /* box comes out too small. We add some code to jiggle things */
     /* in this case. */

{
  int i;
  plc_vector push_up = {{1e-12,1e-12,1e-12}};
  plc_vector push_down = {{.6e-12,0.72e-12,0.69e-12}};

  *bottom = surf->vert_buf[0]; *top = surf->vert_buf[0];

  for(i=0;i<surf->verts;i++) {

    /* Now we check the coordinates, one by one. */

    if (surf->vert_buf[i].c[0] < bottom->c[0]) {

      bottom->c[0] = surf->vert_buf[i].c[0];

    } else if (surf->vert_buf[i].c[0] > top->c[0]) {

      top->c[0] = surf->vert_buf[i].c[0];
      
    }

    if (surf->vert_buf[i].c[1] < bottom->c[1]) {

      bottom->c[1] = surf->vert_buf[i].c[1];

    } else if (surf->vert_buf[i].c[1] > top->c[1]) {

      top->c[1] = surf->vert_buf[i].c[1];
      
    }

    if (surf->vert_buf[i].c[2] < bottom->c[2]) {

      bottom->c[2] = surf->vert_buf[i].c[2];

    } else if (surf->vert_buf[i].c[2] > top->c[2]) {

      top->c[2] = surf->vert_buf[i].c[2];
      
    }

  }

  *top = plc_vect_sum(*top,push_up);
  *bottom = plc_vect_sum(*bottom,push_down);

}


double face_angle_at_vert(surface *surf,int face,int vert)

     /* Procedure computes the angle of the face at vert. */

{
  int ofs,vtf;
  double ang;

  /* First, we check the basics. */

  if (!vert_legal(surf,vert) || !face_legal(surf,face)) {

    fprintf(stderr,"face_angle_at_vert: Bad vert index or face index.\n");
    exit(1);

  }

  ofs = vert_pos_on_face(surf,vert,face);
  
  if (ofs == -1) {

    fprintf(stderr,"face_angle_at_vert: Vert %d is not on face %d!\n",
	    vert,face);

    exit(1);

  }

  /* Now we can work. */

  vtf = surf->vtf[face];
  ofs += vtf;

  plc_vector diffs[2];

  diffs[0] = plc_vect_diff(surf->vert_buf[surf->face_buf[face][(ofs - 1)%vtf]],
			   surf->vert_buf[surf->face_buf[face][(ofs)%vtf]]);
  
  diffs[1] = plc_vect_diff(surf->vert_buf[surf->face_buf[face][(ofs + 1)%vtf]],
			   surf->vert_buf[surf->face_buf[face][(ofs)%vtf]]);

  
  bool ok;
  ang = fabs(plc_angle(diffs[0],diffs[1],&ok));

  if (!ok) { printf("face_angle_at_vert: Can't compute angle, probably this face\n"
		    "                    has repeated vertices."); exit(1); }
  return ang;
}

#ifdef CONVERTED

double total_vertex_angle(surface *surf,int vert)

     /* Procedure computes the total angle at a given vertex. */

{
  int i;
  double total = {0};

  /* First, we check the basics. */
  
  if (!vert_legal(surf,vert)) {

    fprintf(stderr,"total_vertex_angle: Bad vertex index %d.\n",vert);
    exit(1);

  }

  if (surf->incident_angles == NULL) {

    make_incidence_lists(surf);

  }

  /* Now we can work. */

  for(i=0;i<surf->ftv[vert];i++) {

    total += surf->incident_angles[vert][i];

  }

  return total;
}

double gauss_curvature(surface *surf) 

     /* Procedure computes the total Gauss curvature of a surface. */
     /* We get this by computing the angle defect at each vertex.  */

{
  double total = {0};
  int    i;

  for(i=0;i<surf->verts;i++) {

    total += TWO_PI - total_vertex_angle(surf,i);

  }

  return total;

}

#endif  
