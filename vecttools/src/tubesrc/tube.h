
/*

      tube.h : Function prototypes and data types for the standalone
               "tube" application. 

*/

#ifndef __TUBE_H__
#define __TUBE_H__

#if HAVE_CONFIG_H
#include<config.h>
#endif

#ifdef HAVE_ASSERT_H
#include<assert.h>
#endif

#ifdef HAVE_MATH_H
#include<math.h>
#endif

#ifdef HAVE_STDLIB_H
#include<stdlib.h>
#endif

#ifdef HAVE_STRING_H
#include<string.h>
#endif

#include "plCurve.h"
#include "argtable2.h"
#include "uv.h"
#include "../utilib/mangle.h"
#include "../utilib/ordie.h"

#define DEBUG 1 
/* This turns on the various asserts in parts of the code */

int  octrope_error_num;
char octrope_error_str[80];

/**********************************************************************/
/*                                                                    */
/*  Global variables.                                                 */
/*                                                                    */
/**********************************************************************/
 
extern struct arg_dbl *r;
extern struct arg_dbl *g;
extern struct arg_dbl *b;

extern struct arg_dbl *ratio;

extern struct arg_lit *capped;
extern struct arg_lit *closed;
extern struct arg_lit *equalized;
extern struct arg_lit *verbose;


extern struct arg_lit *flared;
extern struct arg_dbl *flare_power;
extern struct arg_dbl *flare_radius;
extern struct arg_dbl *flare_distance;

extern struct arg_lit *bulged;
extern struct arg_dbl *bulge_power;
extern struct arg_dbl *bulge_radius;
extern struct arg_dbl *bulge_start;
extern struct arg_dbl *bulge_end;

extern struct arg_dbl  *radius;
extern struct arg_file *corefile;
extern struct arg_lit  *help;

extern struct arg_end *end;
extern struct arg_end *helpend;

extern plCurve *core;
extern FILE         *infile_fptr,*outfile_fptr;

extern int    nbulges;
extern double *bstarts,*bends,*bpowers,*brads;

extern double DEFAULT_BSTART;
extern double DEFAULT_BEND;
extern double DEFAULT_BPOWER;
extern double DEFAULT_BRAD;

extern int    nflares;
extern double *fpowers, *frads, *fdists;

extern double DEFAULT_FPOWER;
extern double DEFAULT_FRAD;
extern double DEFAULT_FDIST;

extern int MINSIDES;
extern int MAKERINGS;
extern plCurve *rings;

/**********************************************************/
/*                     Surfaces                           */
/**********************************************************/

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

/***************************************************************/
/*                  Function Prototypes                        */
/***************************************************************/ 

void write_surf_to_OFF(surface *surf,FILE *outfile);
int  new_surface(int verts,int faces,surface *surf);

void kill_edge_list(surface *surf);
void kill_incidence_lists(surface *surf);
void kill_surface(surface *surf);
void color_surface(surface *surf,plc_color col);

plc_vector tube_tangent(plCurve *link,int cp,int vt,bool *ok);

plCurve  *plCurve_equalize_sides(plCurve *link,
					   int *target_verts);

plc_vector *random_framezeros(plCurve *core);
void plCurve_bishop_frame(plCurve *link,
			       plc_vector *framezeros,
			       plCurve **frameA, 
			       plCurve **frameB);

void plc_force_closed( plCurve *inLink );
plCurve *split_sharp_corners( plCurve *L );

/* Tube construction code */

void build_ring(int numsteps, double radius, 
		plCurve *L,
		int cp,
		int vt,
		plc_vector frameA, 
		plc_vector frameB,
		surface *tube,
		double s,
		struct uvbuf *uvb);

void match_rings(int startA, int endA, 
		 int startB, int endB,
		 int ab_ofs,
		 plc_color ringcolor,
		 surface *tube); 


int compute_numsteps(plCurve *link, 
		     int cp, int vt, 
		     double rad, 
		     double stepratio);

void close_frame(plCurve *L, int cmp, 
		 plCurve *frameA, 
		 plCurve *frameB, 
		 int numsteps, int *ofs);

surface *make_tube(plCurve *L, 
		   double (*radius)(int comp,int vert),
		   double stepratio,
		   plCurve *frameA,
		   plCurve *frameB,
		   int capped,
		   struct uvbuf *uvb);

#endif
