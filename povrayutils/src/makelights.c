/*******************************************************************

  makelights : This program takes one argument,

                 makelights <file.bbox>

	       and constructs a standard "scientific illustration"
               lighting model for the object in the bounding box. 
               There is some user control on the lighting from the
               command-line.

*******************************************************************/

#include<config.h>

#ifdef HAVE_STDIO_H
#include<stdio.h>
#endif

#ifdef HAVE_STDLIB_H
#include<stdlib.h>
#endif

#ifdef HAVE_MATH_H
#include<math.h>
#endif

#ifdef HAVE_STRING_H
#include<string.h>
#endif

#include "argtable2.h"
#include "plCurve.h"
#include "plsurf.h"

#include "utilib/mangle.h"
#include "bbox.h"
#include "utilib/ordie.h"


/**********************************************************************/

FILE        *infile,*outfile;
char         obj_name[1000];

/**********************************************************************/

void plc_frame(plc_vector x,plc_vector frame[3])

     /* Builds an orthonormal frame with frame[0] = x. */

{
  plc_vector zhat = {{0,0,1}};
  bool ok;
  int i;

  frame[0] = x;
  frame[1] = plc_cross_prod(x,zhat);
  frame[2] = plc_cross_prod(frame[0],frame[1]);

  for(i=0;i<3;i++) {

    frame[i] = plc_normalize_vect(frame[i],&ok);

    if (!ok) { 
      
      fprintf(stderr,"plc_frame: Couldn't build frame on (%g,%g,%g).\n",
	      plc_M_clist(x));

    }

  }

}

void pointlight(FILE *lightfile,plc_vector loc,char *color,bool fill)
{
  fprintf(lightfile,
	  "light_source {\n"
	  " <%g,%g,%g> \n"
	  " color %s\n",
	  plc_M_clist(loc),color);

  if (fill) {

    fprintf(lightfile,"shadowless\n");
    
  }

  fprintf(lightfile,
	  "}\n\n");

}

void cylinderlight(FILE *lightfile,plc_vector location,plc_vector point_at,
		   double radius,char *color,bool fill)

{
  plc_vector frame[3];
  
  /* We first construct a frame for the light. */

  plc_frame(plc_vect_diff(point_at,location),frame);
  frame[1] = plc_scale_vect(radius/2.0,frame[1]);
  frame[2] = plc_scale_vect(radius/2.0,frame[2]);
  
  fprintf(lightfile,
	  "light_source {\n"
	  "<%g,%g,%g>\n"
	  "color %s\n"
	  "cylinder\n"
	  "radius %g\n"
	  "falloff %g\n"
	  "point_at <%g,%g,%g>\n"
	  "area_light <%g,%g,%g>,<%g,%g,%g>,5,5\n"
	  "adaptive 5\n"
	  "jitter\n",
	  plc_M_clist(location),color,radius,2.0*radius,
	  plc_M_clist(point_at),
	  plc_M_clist(frame[1]),plc_M_clist(frame[2]));

  if (fill) {

    fprintf(lightfile,"shadowless\n");

  } 

  fprintf(lightfile,"}\n\n");

}

void kfr_lighting_model(FILE *outfile,double angle,double intensity,plc_vector l,plc_vector u)
{
  plc_vector key,fill;
  plc_vector lookat;
  double dimension;

  lookat = plc_vlincomb(0.5,l,0.5,u);
  dimension = plc_distance(l,u);

  key  = plc_build_vect(10*dimension*cos(angle),-10*dimension*sin(angle),5*dimension);
  fill = plc_build_vect(-key.c[1],key.c[0],key.c[2]);

  fprintf(outfile,
	  "#declare KEYCOLOR = rgb <1,1,1>*%g;\n"
	  "#declare FILLCOLOR = rgb <1,1,1>*%g;\n",
	  intensity,0.5*intensity);

  cylinderlight(outfile,key,lookat,5*dimension,"KEYCOLOR",false);
  cylinderlight(outfile,fill,lookat,5*dimension,"FILLCOLOR",true);
  
}

void bigsoftdome(FILE *outfile,double intensity,plc_vector l,plc_vector u)
{
  surface dome;
  int i;

  dome = icosahedron();
  scale_surface(&dome,5*plc_distance(l,u));
  translate_surface(&dome,plc_vweighted(0.5,l,u));

  fprintf(outfile,
	  "#declare DOMECOLOR = rgb <1,1,1>*%g;\n",
	  intensity);

  for (i=0;i<dome.verts;i++) {

    if (true) {

      pointlight(outfile,dome.vert_buf[i],"DOMECOLOR",true);

    }
 
  }

}

void ambient(FILE *outfile,double intensity)
{
  
  fprintf(outfile,"global_settings { ambient_light rgb<%g, %g, %g> }\n\n",
	  intensity,intensity,intensity);

}

int main(int argc,char *argv[])
{

  double PI = {3.1415926};

  struct arg_dbl *arg_kfrintensity = arg_dbl0("i","intensity","0..1",
						 "strength of main light source");
  struct arg_dbl *arg_domeintensity = arg_dbl0("d","dome","0..1",
						  "strength of \" ambient \" dome light source");
  struct arg_dbl *arg_ambintensity = arg_dbl0("a","ambient","0..1",
						 "strength of ambient light");
  struct arg_file *arg_infile = arg_filen(NULL,NULL,"<file>",1,65000,"input file in .bbox format");
  struct arg_int  *arg_verb = arg_int0("v","verbosity","<0-10>","verbosity");
  struct arg_lit  *arg_help = arg_lit0("?","help","display help message");

  struct arg_end *end = arg_end(20);

  void *argtable[] = {arg_infile,
		      arg_kfrintensity,
		      arg_domeintensity,
		      arg_ambintensity,
		      arg_verb,
		      arg_help,end};
  int nerrors;

  /* Now we parse the arguments */

  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  nerrors = arg_parse(argc,argv,argtable);

  /* special case: '--help' takes precedence over error reporting */
  if (arg_help->count > 0 || argc == 1) {

    printf("makelights constructs POV-ray lights to illuminate the given bounding box.\n");

    fprintf(stderr,"Usage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return 0;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    fprintf(stderr,"\n");
    arg_print_errors(stderr,end,argv[0]);
    fprintf(stderr,"\nUsage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  /* We now start processing the arguments */

  double intensity=0.7,domeintensity,ambintensity = 0.1;

  if (arg_kfrintensity->count > 0) { intensity = arg_kfrintensity->dval[0]; }
  if (arg_domeintensity->count > 0) { domeintensity = arg_domeintensity->dval[0]; } 
  else { domeintensity = (0.15/0.7) * intensity; }
  if (arg_ambintensity->count > 0) { ambintensity = arg_ambintensity->dval[0]; }
  
  int infilenum;

  for(infilenum=0;infilenum < arg_infile->count;infilenum++) {

    plc_vector l,u;

    /* We begin by constructing an outfile */

    outfile = fmangle(arg_infile->basename[infilenum],".bbox",".lights.inc");
    
    /* Now we try to open the input file */

    FILE *infile;

    infile = fopen_or_die(arg_infile->filename[infilenum],"r");
     
    if (!read_bbox(&l,&u,infile)) {
    
      fprintf(stderr,"makelights: %s does not contain valid bbox data.\n",
	      arg_infile->basename[infilenum]);
      exit(1);
      
    }
    
    printf("makelights: Read %g x %g x %g bounding box from file %s.\n",
	   u.c[0] - l.c[0],u.c[1] - l.c[1],u.c[2] - l.c[2],argv[1]);
    
    printf("makelights: Building lighting array.\n");

    kfr_lighting_model(outfile,0.33*PI,intensity,l,u);
    bigsoftdome(outfile,domeintensity,l,u);
    ambient(outfile,ambintensity);

    fclose(outfile);

  }

  exit(0);

}



