#include<config.h>

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_MATH_H
#include <math.h>
#endif

//#include <octrope.h>

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#include <plCurve.h>
#include <argtable2.h>


struct arg_file *corefile;
struct arg_file *output_file;
struct arg_lit  *ring;
struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;
struct arg_end  *end, *helpend;

#ifndef PI
  #define PI 3.1415926
#endif

plCurve *mirror_connect_sum(plCurve *thiscurve,bool add_waist) {


  if (thiscurve == NULL) {

    fprintf(stderr,"mirror_connect_sum: Passed NULL pointer for thiscurve.\n");
    exit(1);

  }

  if (thiscurve->nc > 1) {

    fprintf(stderr,"mirror_connect_sum: Input curve has %d components, this function only works with knots (single component).",thiscurve->nc);
    exit(1);

  }

  if (thiscurve->cp[0].open == true) {

    fprintf(stderr,"mirror_connect_sum: Only works on closed curves.");
    exit(1);

  }

  /* Find a lowest vertex in the Z direction */

  int lowvt,vt;
  double zheight = 1e10;

  for(vt=0;vt<thiscurve->cp[0].nv;vt++) {

    if (thiscurve->cp[0].vt[vt].c[2] < zheight)  {

      lowvt = vt;
      zheight = thiscurve->cp[0].vt[vt].c[2];

    }

  }

  /* Make this the first vertex on the curve. */

  plCurve *working_copy;

  working_copy = plc_copy(thiscurve);

  if (lowvt != 0) {

    for(vt=0;vt<thiscurve->cp[0].nv;vt++) {

      working_copy->cp[0].vt[vt] = thiscurve->cp[0].vt[(vt+lowvt)%thiscurve->cp[0].nv];

    }

    plc_fix_wrap(working_copy);

  }

  /* Now we move the curve rigidly until this first vertex is 3 units above the origin plane */

  plc_vector translate = {{0,0,0}};
  translate.c[0] = -working_copy->cp[0].vt[0].c[0];
  translate.c[1] = -working_copy->cp[0].vt[0].c[1];
  translate.c[2] = 3 - zheight; 
  plc_translate(working_copy,translate);

  /* A good plan is to try to match the resolution of the new sections to the previous */

  int new_verts;
  new_verts = (3.0*plc_num_verts(working_copy))/plc_arclength(working_copy,NULL);
  
  /* We are now ready to copy and start to splice. */

  plCurve *mirroredCurve;
  int nv[2] = {0,0},cc[2] = {0,0}; 
  bool open[2] = {false,false};

  nv[0] = working_copy->cp[0].nv*2 + 4*new_verts;
  nv[1] = 5*(2*PI/3.0)*new_verts;

  mirroredCurve = plc_new(add_waist ? 2:1,nv,open,cc);

  if (mirroredCurve == NULL) {

    fprintf(stderr,"mirror_connect_sum: Couldn't allocate space for %d vertex curve.\n",nv[0]);
    exit(1);

  }

  double zstep = 3.0/(new_verts);
  int i;

  /* Section 1: Climb from plane to curve start. */

  for(i=0,vt=0;i<new_verts;vt++,i++) {

    mirroredCurve->cp[0].vt[vt] = plc_build_vect(working_copy->cp[0].vt[0].c[0],
						 working_copy->cp[0].vt[0].c[1],
						 zstep*i + 0.5*zstep);

  } 

  /* Section 2: Copy the curve */

  for(i=0;i<working_copy->cp[0].nv;vt++,i++) {

    mirroredCurve->cp[0].vt[vt] = working_copy->cp[0].vt[i];

  }

  /* Section 3. Head down to the plane */

  for(i=0;i<new_verts;i++,vt++) {

     mirroredCurve->cp[0].vt[vt] = plc_build_vect(working_copy->cp[0].vt[-1].c[0],
						  working_copy->cp[0].vt[-1].c[1],
						  3 - zstep*i - 0.5*zstep);

  }

  /* Section 4. Continue down to z = -3; */

  for(i=0;i<new_verts;i++,vt++) {

     mirroredCurve->cp[0].vt[vt] = plc_build_vect(working_copy->cp[0].vt[-1].c[0],
						  working_copy->cp[0].vt[-1].c[1],
						  0 - zstep*i - 0.5*zstep);

  }

  /* Section 5: Copy the curve in reverse and flipped over the xy plane */

  for(i=working_copy->cp[0].nv-1;i>=0;vt++,i--) {
    
    mirroredCurve->cp[0].vt[vt] = working_copy->cp[0].vt[i];
    mirroredCurve->cp[0].vt[vt].c[2] *= -1.0;  
    
  }

  /* Section 6: Ascend to the plane */

  for(i=0;i<new_verts;vt++,i++) {
    
    mirroredCurve->cp[0].vt[vt] = plc_build_vect(working_copy->cp[0].vt[0].c[0],
						 working_copy->cp[0].vt[0].c[1],
						 -3 + zstep*i + 0.5*zstep);

  } 

  /* Now add the circle (at radius 5) if we're going to. */

  if (add_waist) {

    double tstep = 2*PI/(mirroredCurve->cp[1].nv);
    double theta;

    for(theta=0,i=0;i<mirroredCurve->cp[1].nv;i++,theta+=tstep) {
			  
      mirroredCurve->cp[1].vt[i] = plc_build_vect(5*cos(theta),5*sin(theta),0);

    }

  }
  
  plc_fix_wrap(mirroredCurve);
 
  /* Clean up and return */

  plc_free(working_copy);

  return mirroredCurve;

}

int main (int argc, char *argv[]) {

  int  nerrors;
  
  void *argtable[] = 
    {
      corefile = arg_filen(NULL,NULL,"<file>",1,65532,"input files"),
      output_file = arg_file0("o","output","<file>","output file"),
      ring = arg_lit0(NULL,"WaistRing","add circle around ``waist'' of mirror knot"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};
  int infilenum;

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("mirrorknot: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"mirrorknot (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"mirrorknot");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"mirrorknot (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      printf("mirrorknot reflects a VECT across the x-y plane to generate symmetric composites.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  /* Now we iterate over input files */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    plCurve *thiscurve;
    FILE *infile_fptr;
    
    /* We begin by loading the current core link. */
    
    infile_fptr = fopen(corefile->filename[infilenum],"r");
    
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"mirrorknot: Couldn't open file %s.\n",corefile->filename[infilenum]);
      continue;  /* Try the next file */
      
    }
  
    int plr_error_num;
    char plr_error_str[1024];

    thiscurve = plc_read(infile_fptr,
			 &plr_error_num,plr_error_str,sizeof(plr_error_str));
    
    if (plr_error_num > 0) {   /* This is the signal for an error. */
    
      fprintf(stderr,"mirrorknot: link reading error\n%s\n",plr_error_str);
      continue;  /* Try the next file */
      
    }

    /* We have now succeeded in loading a file */

    plCurve *mirroredCurve;

    mirroredCurve = mirror_connect_sum(thiscurve,ring->count > 0);

    /* Now we construct a new filename and save it */

    FILE *outfile;
    char outfile_name[1000];

    if (infilenum < output_file->count) { 

      outfile = fopen(output_file->filename[infilenum],"w");     
      strncpy(outfile_name,output_file->filename[infilenum],sizeof(outfile_name));

    } else {

      char basename[512];
      char *endofname;
      char mirrorname[1000];
      
      endofname = strcasestr(corefile->basename[infilenum],".vect");
      strncpy(basename,corefile->basename[infilenum],endofname - corefile->basename[infilenum]);
      sprintf(mirrorname,"%s_#_%sm.vect",basename,basename);
      
      outfile = fopen(mirrorname,"w");
      strncpy(outfile_name,mirrorname,sizeof(outfile_name));
      
    }
    
    if (outfile == NULL) {
      
      fprintf(stderr,"mirrorknot: Couldn't open output file %s.\n",outfile_name);
      exit(1);

    } 

    plc_write(outfile,mirroredCurve);

    printf("mirrorknot: Wrote %d component, %d vertex mirrored knot to %s.\n",
	   mirroredCurve->nc,plc_num_verts(mirroredCurve),outfile_name);

    fclose(outfile);

    plc_free(thiscurve);
    plc_free(mirroredCurve);

  }

  exit(0);
    
}


