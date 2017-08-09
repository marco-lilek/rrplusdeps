/*

    plsurf.h : The definitions and declarations for surface objects.
               A surface object, defined from Geomview's OFF format,
	       is a collection of polygons, possibly with shared vertices.

*/

#ifndef __SURFACES_H__
#define __SURFACES_H__

#include "plCurve.h"

extern const plc_color SURFACE_RED;
extern const plc_color SURFACE_GREEN;
extern const plc_color SURFACE_BLUE;
extern const plc_color SURFACE_GRAY;

#define NOT_COMPUTED -1
#define NO_FACE -1

typedef struct {

  int verts;
  int faces;

  plc_vector *vert_buf;  /* All the vertices, numbered 0 .. verts-1. */
  int            *vtf;       /* # of vertices on the Nth face. */
  int           **face_buf;  /* A list of the vertices on the nth face. */
  plc_color  *col_buf;   /* A list of the colors of each face. */

  int             edges;  
  int           **edge_list;    

  /* This data structure gives a list of the edges in the surface in  */
  /* the following format:                                            */
  
  /*  <#faces edge appears in> <v_1> <v_2> <f_1> <o_1> ...            */

  /* where v_1 and v_2 are the start and end of the edge, the f_i     */
  /* are face numbers, and the <o_i> are edge orientations.           */
  
  /* This is a computed data structure, filled by calling make_edge_list */

  int     *ftv;

  /* This is an array of the number of faces adjacent to each vertex. */

  int     **incident_faces;

  /* This data structure gives a list of the faces incident to each */
  /* vertex of the surface, in the following form. */

  /* <f_1> <f_2> ... */

  /* where the f_i are face numbers. These appear in a particular order:
     the faces go in _counterclockwise_ order around the vertex. If the 
     vertex is a boundary vertex, special rules apply:

            1) The incident faces are listed in counterclockwise order
               from the first to the last. 

            2) An extra face, denoted NO_FACE (= -1) is added to the 
               end of the list for this vertex. 

            3) An extra angle, 0, is added to incident_angles.

            4) ftv is incremented by 1, to take into account the extra
               face.  

     All of these modifications are picked up and handled by the other
     procedures (namely, tvector_angle, and tvector_at_angle) which need
     them. */

  /* This is a computed data structure, filled by make_incidence_lists. */
  /* The number of entries in the array is given by ftv. */

  double  **incident_angles;

  /* This data structure, used by exp, contains the angle at the corner */
  /* of the face <f_i> at the vertex v_j, where the face indexing comes */
  /* from the incident_faces list above. Data format is: */

  /* <a_1> <a_2> ... */

  /* This is also filled by make_incidence_lists. */
  /* The number of entries is given by ftv. */

} surface;

struct uvbuf {

  int verts;
  plc_vector *uv;  /* uv[i].c[0] = u, uv[i].c[1] = v, uv[i].c[2] = 0 */ 

};

/**********************************************/
/*               Surface I/O                  */
/**********************************************/

int skip_whitespace_and_comments(FILE *infile);

int scandoubles(FILE *infile,int ndoubles, ... );
int scanints(FILE *infile,int nints, ... );
int scancolor(FILE *infile,plc_color *thiscolor);

int  load_surf_from_OFF(surface *surf,FILE *infile);
void write_surf_to_OFF(surface *surf,FILE *outfile);
void write_surf_to_POV(surface *surf,FILE *outfile,struct uvbuf *uvb);

void print_edge_list(surface *surf);
void print_incidence_lists(surface *surf);
void print_shared_faces(surface *surf);

/**********************************************/
/*          u-v texture coordinates           */
/**********************************************/

void write_uvfile(struct uvbuf *uvb,FILE *outfile);
int  load_uvfile(struct uvbuf *uvb,FILE *infile);

/**********************************************/
/*            Surface Maintenance             */
/**********************************************/

int     new_surface(int verts,int faces,plc_color def_col,surface *surf);
void    new_face(surface *surf,int face_num,int vtf,...);

void    kill_surface(surface *surf);
void    kill_edge_list(surface *surf);
void    kill_incidence_lists(surface *surf);

void    reindex_face(surface *surf,int face,int start_vert);

//void    triangulate_surface(surface *surf);
//void    subdivide_edges(surface *surf,int num_edges,int *edges);
void    subdivide_faces(surface *surf);

//surface consolidate_surface(surface *surf,double precision);
//void    orient_surface(surface *surf);

/*********************************************/
/*        Surface Operations                 */
/*********************************************/

void    scale_surface(surface *surf,double scale);
//void    normalize_surface_volume(surface *surf);
void    translate_surface(surface *surf,plc_vector vec);
void    spherical_projection_surface(surface *surf);

//double  vert_cone_angle(int vert,surface *cone);
//void    unroll_cone(surface *cone);

//surface join_surfaces(int n_surfs, ... );
//surface pushoff_surface(surface *surf,double eps);

//void     reverse_surface(surface *surf);
//surface *unfold_surface(surface *surf,int face,int edge);

//void     close_surface_by_cone(surface *surf,vector vertex);

#ifdef CONVERTED

/**********************************************/
/*        Surface Components                  */
/**********************************************/

int      is_edge_on_bdy(surface *surf,int edge);
int      is_vert_on_bdy(surface *surf,int vert);

int      get_boundary_edges(surface *surf,int **boundary_edges);
int      get_boundary_verts(surface *surf,int **boundary_verts);

#endif

/**********************************************/
/*         Surface Integrity Checks           */
/**********************************************/

//int  sanity_surface(surface *surf);

bool  is_triangulated(surface *surf);

//int  is_top_mfld(surface *surf);
//int  is_oriented(surface *surf);
//int  is_closed(surface *surf);

//int  edge_list_ok(surface *surf);
//int  verts_unique(surface *surf);
bool  faces_ok(surface *surf);
//int  faces_planar(surface *surf);
//int  colors_good(surface *surf);

bool  face_ok(surface *surf,int i);
//int  face_planar(surface *surf,int i);
bool  goodcolor(plc_color col);

bool  vert_legal(surface *surf,int vert);
bool  edge_legal(surface *surf,int edge);
bool  face_legal(surface *surf,int face);


/***********************************************/
/*       Computing Surface Geometry            */
/***********************************************/

//int    inside(surface *surf,vector *point);

//double gauss_surface(surface *surf,vector *origin);
//double polygon_proj(vector *origin, int num_sides, vector *verts);
//double triangle_proj(vector *origin, vector *a, vector *b, vector *c);

//double area_surface(surface *surf);
//double volume_surface(surface *surf);
//double area_face(surface *surf,int i);
//double edge_length(surface *surf,int edge);

//void   face_outline(surface *surf,int face,polyline *pline);

plc_vector face_normal(surface *surf,int i);
//vector edge_normal(surface *surf,int i);
plc_vector vert_normal(surface *surf,int i);

//vector inward_normal(surface *surf,int edge,int face);
//vector outward_normal(surface *surf,int edge,int face);

//vector barycenter(surface *surf,int face);
plc_vector vertex_center_of_mass(surface *surf);
void   bbox_surface(surface *surf,plc_vector *bottom,plc_vector *top);

double face_angle_at_vert(surface *surf,int face,int vert);
//double total_vertex_angle(surface *surf,int vert);
//double gauss_curvature(surface *surf);

#ifdef CONVERTED

/**********************************************/
/*      Implicit Surface Construction         */
/**********************************************/

#define FIND_ZERO   1
#define FIND_LINEAR 2

vector find_zero(vector A, vector B,int iterations,
		 double (*function)(vector loc,void *external_data),
		 void *external_data);

vector find_linear(vector A, vector B,double(*function)(vector loc,void *external_data),
		   void *external_data);

surface implicit_surface(double (*function)(vector loc,void *external_data),
		  	 vector bbox[2],int gridsteps,void *external_data,int mode);

#endif

/**********************************************/
/*      Standard Surfaces                     */
/**********************************************/

surface icosahedron();
surface geodesic_dome(int subdivision,int recursion);

/***********************************************/
/*      Edge and Incidence List Code           */
/***********************************************/

struct face_node {

  struct face_node *next;
  int    face;
  int    orientation;

};

struct elist_node {

  struct elist_node *next;
  int                vertex;
  struct face_node  *root;
  int                num_faces;

};

int  edges_surface(surface *surf);
void make_edge_list(surface *surf);
void make_incidence_lists(surface *surf);

int  face_pos_in_elist_record(surface *surf,int face,int edge);

/*********************************************/
/*     Structure of the Triangulation        */
/*********************************************/

int  edge_shared_faces(surface *surf,int edge_one,int edge_two,int **result);
int  vert_shared_faces(surface *surf,int vert_one,int vert_two, int **result);
int  edge_vert_shared_faces(surface *surf,int edge,int vert,int **result);

bool  vert_on_face(surface *surf,int vert,int face);
bool  edge_on_face(surface *surf,int edge,int face);

int  vert_pos_on_face(surface *surf,int vert,int face);
int  least_vert_pos_on_face(surface *surf,int face);
int  face_pos_at_vert(surface *surf,int face,int vert);
int  face_pos_in_edge_rec(surface *surf,int face,int edge);
int  next_vert_on_face(surface *surf,int vert,int face);

int  vert_pos_opposite_edge_on_face(surface *surf,int edge,int face);
void face_edges_meeting_at_vert(surface *surf,int face,int vert,int *edges);

bool  faces_adjacent(surface *surf,int face_one,int face_two,int *edge);
int  faces_vertex_adjacent(surface *surf,int f_one,int f_two,int **verts);
void adjacent_faces_to_face(surface *surf,int face,int **faces); 

int  edge_number(surface *surf,int vert_one,int vert_two);
int  adjacent_face(surface *surf,int edge,int face);

void edge_positions_on_face(surface *surf,int edge,int face,int edgepos[2]);
int  edge_orientation_on_face(surface *surf,int edge,int face);
void edges_on_face(surface *surf,int face,int **edges);
int  next_edge_on_face(surface *surf,int face,int edge);

void edges_at_vert(surface *surf,int vert,int *edges);

/************************************************/
/*             Surface Color Models             */
/************************************************/

void color_surface(surface *surf,plc_color col);
//void color_by_proximity(surface *surf,polyline *curve,double radius);

#ifdef CONVERTED 

/*************************************************/
/*                 Surfpt Code                   */
/*************************************************/

/* The next data structure represents a point on a TRIANGULATED surface,
   given in barycentric coordinates by face number. There are three types
   of surfpt's:

   FACEPT
   EDGEPT
   VERTPT

   with corresponding data. We require that the edge list be computed,
   before we can start placing surfpts on a surface. 

   By barycentric coordinates, we mean that the point is located at:

        floc[0] * vert[0] + floc[1] * vert[1] + floc[2] * vert[2],

   and that the sum of the floc[i] is equal to one. Notice that the 
   meaning of these coordinates depends explicitly on the ordering of
   the faces in the triangle!

   When we define the position of an EDGEPT, we mean that it is at:

        eloc * surf->edge_list[edge][1] + (1-eloc) surf->edge_list[edge][2]

*/

enum surfpt_type { FACEPT, EDGEPT, VERTPT, BADPT };

typedef struct {

  surface     *surf;
  
  enum surfpt_type type;

  int         face;
  double      floc[3];

  int         edge;
  double      eloc;

  int         vert;

} surfpt;

int     surfpt_legal(surfpt p);

surfpt  tofacept(surfpt p,int face);
surfpt  toedgept(surfpt p,int edge);
surfpt  tovedgept(surfpt p);

vector  tovector(surfpt p);
int     fromvector(surface *surf,int face,vector v,surfpt *result);

char    *print_surfpt(surfpt p);
surfpt  read_surfpt(char *s);

double  surfpt_total_angle(surfpt p);
double  surfpt_distance(surfpt a,surfpt b);

surfpt  surfpt_interpolate(surfpt a,surfpt b,double t);

int     is_on_face(surfpt a,int face);
int     on_adjacent_faces(surfpt a,surfpt b,
			  int *face_a,int *face_b,
			  int *edge);

int     on_vertex_adjacent_faces(surfpt a,surfpt b,
				 int **a_faces,int **b_faces,
				 int **shared_verts_this_pair,
				 int **shared_verts);

int     at_same_vertex(surfpt a,surfpt b,int *vert);
int     on_same_edge(surfpt a,surfpt b,int *edge); 
int     on_same_face(surfpt a,surfpt b,int *face);
  
int     on_faces(surfpt a,int **result);
int     on_bdy(surfpt a);
int     same_surfpt(surfpt a,surfpt b);

int     pair_on_faces(surfpt a,surfpt b,int **result);

/******************************************************/
/*                Tvector Code                        */
/******************************************************/

/* The last data structure represents a direction in the tangent space */
/* to the surface. It is written in barycentric coordinates on the face, */
/* i.e. the tangent direction is given by 

   floc[0] * vert[0] + floc[1] * vert[1] + floc[2] * vert[2],

   and the sum of the floc[i] is equal to ZERO. All tvectors are valid
   at a surfpt of type FACEPT. At a point of type EDGEPT or VERTPT, special
   rules apply:

         1) We think of the tangent space at an EDGEPT as divided into 
            two halves, each containing 180 degrees worth of angle, pointing
            _in_ to the face in question.
	    
	    Therefore, at an EDGEPT, it makes sense to speak of some tvectors
	    as lying on each face... procedures like tvector_at_angle and 
	    change_chart are capable of switching faces as required.

	 2) We think of the tangent space at a VERTPT as divided into sections,
	    each one along the corner of the faces which meet at the vertex,
	    containing the directions pointing _in_ to the face in question.

	    We notice that in general, there are not a total of PI radians
	    of angle at a vertpt. We deal with this situation using the 
	    procedure surfpt_total_angle, which returns the angle sum at
	    a given vertpt. 

	    It will sometimes be convenient to discuss absolute angle measures
	    (such as PI/2) at these points; at other times, we will measure
	    "halfway around" the angle. 

	    Procedures like "change_chart" are quite capable of switching 
	    faces on a tvector as required.
*/

typedef struct {

  double floc[3];
  int    face;

  surface *surf;

} tvector;

int     tvector_legal(tvector v);

vector  tvector_tovector(tvector v);
tvector tvector_fromvector(surface *surf,int face,vector v); 

tvector tvector_from_vedgept(surfpt p,int face);
tvector tvector_to_vedgept(surfpt p,int face);
tvector tvector_inward_normal_at_vedgept(surfpt p,int face);
tvector tvector_outward_normal_at_vedgept(surfpt p,int face);

double  tvector_angle(tvector v,tvector w,surfpt p);
tvector tvector_at_angle(surfpt p,tvector t,double theta);

double  tvector_norm(tvector v);
void    tvector_normalize(tvector *t);
tvector tvector_scalarmultiply(double k,tvector *t);

tvector tvector_perp(surfpt p,tvector v);
matrix  face_basis_from_tvector(tvector v_one);

tvector tvector_interpolate(surfpt loc,tvector a,tvector b,double t);

/*****************************************************/
/*                  Surfpline Code                   */
/*****************************************************/

/* It will sometimes be convenient to build a linked list of surfpts. */

struct surfpline_node {

  surfpt pt;
  struct surfpline_node *next;

};

/* But a general surface polyline uses the usual buffer. */

typedef struct {

  surface *surf;
  
  int verts;
  int closed;       /* If this tag is set to TRUE, then we check that 
		       the last vertex can be connected to the first
		       as a valid surfpline. */

  surfpt *vert_buf;
  
} surfpline;

int     surfpline_legal(surfpline *surfpl);

void    llist_to_surfpline(struct surfpline_node *root,surfpline *surfpl);
void    kill_llist(struct surfpline_node **root);

void    new_surfpline(int verts,surfpline *surfpl);
void    kill_surfpline(surfpline *surfpl);

void    topolyline(surfpline *surfpl,polyline *result);

tvector tangent_vector_surfpline(surfpline *surfpl,int vert);

void    reverse_surfpline(surfpline *surfpl);
void    close_surfpline(surfpline *surfpl);

surfpline *join_surfplines(int num_plines,surfpline **join_list);
surfpline *surfpline_pushoff(surfpline *surfpl,double length);
surfpline *surfpline_ribbon(surfpline *surfpl,double width);

/*************************************************/
/*       Slicing Surfaces along Surfplines       */
/*************************************************/

typedef struct crossing_rec_struct {

  int    vert;
  double loc;

  struct crossing_rec_struct *next;

} crossing_rec;
  
typedef struct edge_rec_struct {

  int    ends[2];
  int    original;

  crossing_rec           *cr_root;
  struct edge_rec_struct *next;

} edge_rec;

void list_segments_on_face(surface *surf,surfpline *surfpl,
			   edge_rec **er_buf,surfpt *new_vert_buf,
			   int *verts);

void find_crossings_on_face(edge_rec *this_face,
			    surfpt *new_vert_buf,
			    int *verts);

void assemble_new_faces_on_face(edge_rec *this_face,
				int *new_face_buf,
				int *faces);

void slice_surface(surface *surf,surfpline *surfpl,surface *new_surf);

/*************************************************/
/*            Geodesics on Surfaces              */
/*************************************************/

double  face_length(surfpt a,tvector t);

tvector change_chart(surfpt p,tvector t);

void    surface_exp(surfpt start, tvector dir,double length,
		    surfpline *geodesic);

surfpt  surface_exp_target(surfpt start,tvector dir,double length);

int     short_geodesic(surfpt start,surfpt end,struct surfpline_node **geo);
int     triangle_intersections_with_edges(surfpt  vertex_pt,
					  tvector tv_start,
					  tvector tv_end,
					  double  theta,
					  struct surfpline_node **root);

#endif
#endif


  






