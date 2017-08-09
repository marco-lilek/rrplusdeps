/*******************************************************************

  makecam : This program takes one argument,

              makecam <file.bbox>

	    and constructs a POV-ray camera aimed at the box 
            described in the file. There are three camera modes:
            print, screen, and video, each choosing camera data
            to look right at different viewing distances.

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
#include "utilib/mangle.h"
#include "bbox.h"
#include "utilib/ordie.h"

double ASPECTRATIO = (4.0/3.0);


/**********************************************************************/

FILE        *infile,*outfile;
char         obj_name[1000];

/**********************************************************************/

void anglecam(FILE *camfile,double angle,double buffer,plc_vector l,plc_vector u)

     /* Creates a perspective camera with a given angle of view which
	displays the entire box (l,u) with a buffer expressed as a
	fraction of the largest dimension of the box. 

	We view the box so that the y-z plane is projected on the screen. */

{
  plc_vector lookat;
  plc_vector loc;

  lookat = plc_vweighted(0.5,l,u);
  
  double dist,distz,disty,DtoR = (3.1415926/180);

  distz = ASPECTRATIO*0.5*(u.c[2]-l.c[2])*(1+buffer)/tan(angle*DtoR/2.0) + (u.c[1] - l.c[1])/2.0;
  disty = 0.5*(u.c[1]-l.c[1])*(1+buffer)/tan(angle*DtoR/2.0) + (u.c[2] - l.c[2])/2.0;
  dist = (distz > disty) ? distz : disty;

  plc_vector dir = {{1,0,0}};

  loc = plc_vlincomb(1,lookat,dist,dir);

  fprintf(camfile,
	  "camera { \n"
	  "  perspective \n"
	  "  angle    %g \n"
	  "  rotate   <90,0,0> \n"
	  "  location <%g,%g,%g> \n"
	  "  look_at  <%g,%g,%g> \n"
	  "}\n",
	  angle,
	  loc.c[0],loc.c[2],-loc.c[1],
	  lookat.c[0],lookat.c[2],-lookat.c[1]);

  /* Since we are rotating the coordinate system, we need to rotate
     the location and look_at vectors as well. */

}
	   
void screen_camera(FILE *outfile,plc_vector l,plc_vector u)
{
  anglecam(outfile,20,0.2,l,u);
}
  
void print_camera(FILE *outfile,plc_vector l,plc_vector u)
{
  anglecam(outfile,10,0.1,l,u);
}

void video_camera(FILE *outfile,plc_vector l,plc_vector u)
{
  anglecam(outfile,26,0.1,l,u);
}

int main(int argc,char *argv[])
{
  struct arg_file *arg_infile = arg_filen(NULL,NULL,"<file>",1,65000,"input file in .bbox format");
  struct arg_dbl  *arg_aspect = arg_dbl0("a","aspectratio","<real>","horiz to vert aspect ratio (default 4/3)");
  struct arg_int  *arg_verb = arg_int0("v","verbosity","<0-10>","verbosity");
  struct arg_lit  *arg_help = arg_lit0("?","help","display help message");

  struct arg_end *end = arg_end(20);

  void *argtable[] = {arg_infile,
		      arg_aspect,
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

    printf("makecam constructs a POV-ray camera aimed at the bounding box\n"
	   "given in <file.bbox>. There are options for print, screen and video.\n");

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

  if (arg_aspect->count > 0) { 
    
    ASPECTRATIO = arg_aspect->dval[0];

  }

  int infilenum;

  for(infilenum=0;infilenum < arg_infile->count;infilenum++) {

    plc_vector l,u;

    /* We begin by constructing an outfile */

    outfile = fmangle(arg_infile->basename[infilenum],".bbox",".camera.inc");
    
    /* Now we try to open the input file */

    FILE *infile;

    infile = fopen_or_die(arg_infile->filename[infilenum],"r");
     
    if (!read_bbox(&l,&u,infile)) {
    
      fprintf(stderr,"makecam: %s does not contain valid bbox data.\n",
	      arg_infile->basename[infilenum]);
      exit(1);
      
    }
    
    printf("makecam: Read %g x %g x %g bounding box from file %s.\n",
	   u.c[0] - l.c[0],u.c[1] - l.c[1],u.c[2] - l.c[2],argv[1]);
    
    printf("makecam: Building print, screen, and video cameras.\n");

    fprintf(outfile,"#declare print_camera = \n");
    print_camera(outfile,l,u);

    fprintf(outfile,"#declare screen_camera = \n");
    screen_camera(outfile,l,u);

    fprintf(outfile,"#declare video_camera = \n");
    video_camera(outfile,l,u);
    
    fclose(outfile);

  }

  exit(0);

}



