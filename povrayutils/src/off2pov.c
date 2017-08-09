/*******************************************************************

  off2pov : This program takes one argument,

              off2pov <file>

	      and converts surface geometry in the OFF format to 
	      a POVray smooth triangle mesh. We expect the input 
	      filename to be in the form <name.off>, and write the 
	      output filename in the form <name.inc>, for inclusion
	      in a .pov scene description file. We don't include 
	      a camera, or any lights, so the resulting file is 
	      not "ready-to-render" in any sense. 

	      The object is given the name 'name' in the inc file.

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

#ifdef HAVE_ASSERT_H
#include<assert.h>
#endif

#include "argtable2.h"
#include "plCurve.h"
#include "plsurf.h"

#include "utilib/mangle.h"
#include "bbox.h"
#include "utilib/ordie.h"

/**********************************************************************/

surface      surf;
FILE        *infile,*outfile,*bboxfile;
char         obj_name[1000];

/**********************************************************************/

int main(int argc,char *argv[])
{
  plc_vector bbox[2];

  struct arg_file *arg_infile = arg_filen(NULL,NULL,"<file.off>",1,65000,"input file in OFF format");
  struct arg_file *arg_outfile = arg_filen("o","output","<file>",0,65000,"output filename");
  struct arg_lit  *arg_uvfile = arg_lit0(NULL,"uvfile","include uv coordinates for "
					 "texture mapping from file");
  struct arg_file *arg_uvfilename = arg_file0(NULL,"uvfilename","<file.uv>","uv texture coordinate file");

  struct arg_int  *arg_verb = arg_int0("v","verbosity","<0-10>","verbosity");
  struct arg_lit  *arg_help = arg_lit0("?","help","display help message");

  struct arg_end *end = arg_end(20);

  void *argtable[] = {arg_infile,
		      arg_outfile,
		      arg_uvfile,
		      arg_uvfilename,
		      arg_verb,arg_help,end};
  int nerrors;
  struct uvbuf uvb;

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

    printf("off2pov converts a Geomview OFF file into a POV-ray mesh\n"
	   "suitable for inclusion as a .inc file into a POV-ray scene.\n"
	   "No camera or lights are provided, so the resulting scene is\n"
	   "NOT ready to render in any sense.\n");

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

  int infilenum;

  for(infilenum=0;infilenum < arg_infile->count;infilenum++) {

    /* We begin by constructing an outfile filename */

    char outfile_name[1024];
    
    if (arg_outfile->count == 0) {
      
      nmangle(outfile_name,sizeof(outfile_name),arg_infile->basename[infilenum],".off",".inc");
      bboxfile = fmangle(arg_infile->basename[infilenum],".off",".bbox");

    } else {

      strncpy(outfile_name,
	      arg_outfile->filename[infilenum % arg_outfile->count],
	      sizeof(outfile_name));
      char *bboxname;

      bboxname = mangle(outfile_name,".inc",".bbox");
      bboxfile = fopen_or_die(bboxname,"w");
      free(bboxname);
      
    }

    /* Now we try to open the outfile. */
     
    FILE *outfile;

    outfile = fopen_or_die(outfile_name,"w");
     
    /* Now we try to open the input file */

    FILE *infile;

    infile = fopen_or_die(arg_infile->filename[infilenum],"r");
     
    /* Now we try to load the OFF surface from the input file. */
    
    if (!load_surf_from_OFF(&surf,infile)) {
      
      fprintf(stderr,"off2pov: %s does not contain valid OFF data.\n",
	      arg_infile->basename[infilenum]);
      exit(1);
      
    }
    
    printf("off2pov: Read %d vert, %d face surface from file %s.\n",
	   surf.verts,surf.faces,argv[1]);

    if (arg_uvfile->count > 0) {

      FILE *uvfile;
      char *uvname;

      if (arg_uvfilename->count > 0) {
	
	uvfile = fopen_or_die(arg_uvfilename->filename[0],"r");

      } else {

	uvname = mangle(arg_infile->filename[infilenum],".off",".uv");
	uvfile = fopen_or_die(uvname,"r");
	free(uvname);

      }
	
      if (!load_uvfile(&uvb,uvfile)) {

	printf("off2pov: Couldn't load uvfile.\n");
	exit(1);

      }

      assert(uvb.verts == surf.verts);

    }

    fprintf(outfile,"/* Object '%s' */\n",arg_infile->basename[infilenum]);
    
    if (arg_uvfile->count > 0) {
      write_surf_to_POV(&surf,outfile,&uvb);
    } else {
      write_surf_to_POV(&surf,outfile,NULL);
    }
    fclose(outfile);
    
    printf("off2pov: Smooth triangle mesh written to %s\n",
	   outfile_name);

    bbox_surface(&surf,&bbox[0],&bbox[1]);
    printf("off2pov: Surface bounding box: (%g,%g,%g,)-",plc_M_clist(bbox[0]));
    printf("(%g,%g,%g).\n",plc_M_clist(bbox[1]));
    write_bbox(bbox[0],bbox[1],bboxfile);
    fclose(bboxfile);
    
    kill_surface(&surf);


  }

  exit(0);

}



