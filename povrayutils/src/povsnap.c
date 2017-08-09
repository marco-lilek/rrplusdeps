/*******************************************************************

  povsnap : This program takes one argument,

                 povsnap <file.off>

	    converts it to povray format, builds camera and lights, 
            and renders it to a PNG file, cropping the PNG file using
            ImageMagick.

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

#ifdef HAVE_SYS_STAT_H
#include<sys/stat.h>
#endif

#include "argtable2.h"
#include "plCurve.h"
#include "plsurf.h"

#include "utilib/mangle.h"
#include "bbox.h"
#include "utilib/ordie.h"

int QUIET = 0;

/**********************************************************************/

FILE        *infile,*outfile;
char         obj_name[1000];

/**********************************************************************/

int main(int argc,char *argv[])
{

  struct arg_lit  *arg_print = arg_lit0("p","print","generate at print resolution (1024x768)");
  struct arg_lit  *arg_screen = arg_lit0("s","screen","generate at screen resolution (320x200)");
  struct arg_lit  *arg_video = arg_lit0("v","video","generate at video resolution (640x480)");
  struct arg_lit  *arg_texture = arg_lit0("t","texture","look for texture data");
  struct arg_lit  *arg_display = arg_lit0(NULL,"NoDisplay","turn off povray display");
  struct arg_int  *arg_axis = arg_int0("a","axis","<1-3>","orient to axis of inertia");
  struct arg_lit  *arg_science = arg_lit0(NULL,"scienceview","orient to default 'scientific' point of view");
  struct arg_dbl  *arg_spin = arg_dbln(NULL,"spin","<rad>",0,3,"spin around each axis (may be repeated, only works for vects)");
  struct arg_dbl  *arg_tuberadius = arg_dbl0(NULL,"TubeRadius","<radius>","radius for tube around VECT files");

  struct arg_file *arg_infile = arg_filen(NULL,NULL,"<file>",1,65000,"input file in .bbox format");
  struct arg_int  *arg_verb = arg_int0("v","verbosity","<0-10>","verbosity");
  struct arg_lit  *arg_help = arg_lit0("?","help","display help message");
  struct arg_lit  *arg_quiet = arg_lit0("q","quiet","run in batch mode");

  struct arg_end *end = arg_end(20);

  void *argtable[] = {arg_infile,
		      arg_quiet,
		      arg_print,
		      arg_screen,
		      arg_video,
		      arg_axis,
		      arg_science,
		      arg_spin,
		      arg_tuberadius,
		      arg_texture,
		      arg_display,
		      arg_verb,
		      arg_help,end};
  int nerrors;
  double tube_radius = 0.5;

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

    printf("povsnap creates a snapshot of an object given as a Geomview OFF or VECT using povray.\n");

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

  if (arg_quiet->count > 0) {

    QUIET = 1;

  }

  if (arg_tuberadius->count > 0) {

    tube_radius = arg_tuberadius->dval[0];

  }

  int infilenum;

  for(infilenum=0;infilenum < arg_infile->count;infilenum++) {

    char basename[2048];

    if (strstr(arg_infile->basename[infilenum],".vect") != NULL) { /* This is a vect! */

      if (arg_axis->count > 0) { /* Need to call orientvect on this guy */

	printf("povsnap: Orienting vect file to axis %d...",arg_axis->ival[0]);

	char cmdline[1024];

	if (arg_spin->count == 0) {
	  sprintf(cmdline,"orientvect -a %d -o temp.vect %s",arg_axis->ival[0],arg_infile->basename[infilenum]);
	} else {
	  double spin[3] = {0,0,0};
	  int i;
	  for(i=0;i<arg_spin->count;i++) { spin[i] = arg_spin->dval[i]; }

	  sprintf(cmdline,"orientvect -a %d -o temp.vect -s %g -s %g -s %g %s",arg_axis->ival[0],spin[0],spin[1],spin[2],arg_infile->basename[infilenum]);
	}

	system(cmdline);
	printf("done\n");

	printf("povsnap: Tubing oriented vect file at radius %g...",tube_radius);
	sprintf(cmdline,"tube -r %g temp.vect; rm temp.vect",tube_radius);
	system(cmdline);
      
	char *temp;
	temp = mangle(arg_infile->basename[infilenum],".vect",".tube.off");
	strcpy(basename,temp);
	free(temp);
	
	sprintf(cmdline,"mv temp.tube.off %s",basename);
	system(cmdline);

      } else if (arg_science->count > 0) {  

	printf("povsnap: Orienting vect file to default 'scientific' view...");

	char cmdline[1024];

	sprintf(cmdline,"orientvect -a 1 --matrix=\"{{0.77791, 0.583432, 0.233373},{-0.6, 0.8, 0.},{-0.186698,-0.140024, 0.972387}}\" "
		"-o temp.vect  %s",arg_infile->basename[infilenum]);	
	system(cmdline);

	printf("done\n");

	printf("povsnap: Tubing oriented vect file at radius %g...",tube_radius);
	sprintf(cmdline,"tube -r %g temp.vect; rm temp.vect",tube_radius);
	system(cmdline);
      
	char *temp;
	temp = mangle(arg_infile->basename[infilenum],".vect",".tube.off");
	strcpy(basename,temp);
	free(temp);
	
	sprintf(cmdline,"mv temp.tube.off %s",basename);
	system(cmdline);

      } else {

	printf("povsnap: Tubing vect file at radius %g...",tube_radius);
	char cmdline[1024];
	sprintf(cmdline,"tube -r %g %s",tube_radius,arg_infile->basename[infilenum]);
	system(cmdline);
	
	char *temp;
	temp = mangle(arg_infile->basename[infilenum],".vect",".tube.off");
	strcpy(basename,temp);
	free(temp);

	printf("done\n");

      }

    } else {

      strcpy(basename,arg_infile->basename[infilenum]);

    }	

    char *camname,*lightname,*objname,*bboxname,*dirname,*scenename, *pngname;

    /* We begin by constructing output names. */

    dirname   = mangle(basename,".off","");
    
    objname   = mangle(basename,".off",".inc");
    bboxname  = mangle(objname,".inc",".bbox");
    camname   = mangle(bboxname,".bbox",".camera.inc");
    lightname = mangle(bboxname,".bbox",".lights.inc");
    
    scenename = mangle(basename,".off",".pov");
    pngname = mangle(scenename,".pov",".png");

    /* Now we create a directory for our output */
    
    char cmdstring[1024];
    sprintf(cmdstring,"rm -fr %s",dirname);
    system(cmdstring);
    mkdir_or_die(dirname,S_IRWXU | S_IRGRP | S_IXGRP);
    chdir_or_die(dirname);
    
    /* We first create a scene file. */

    if (!QUIET) { printf("povsnap: Building scene file..."); }

    FILE *scenefile;

    scenefile = fopen_or_die(scenename,"w");

    fprintf(scenefile,
	    "#include \"finish.inc\"\n"
	    "#include \"colors.inc\"\n"
	    "#include \"textures.inc\"\n" 
	    "#include \"%s\" \n"
	    "#include \"%s\" \n\n",
	    camname,lightname);

    fprintf(scenefile,
	    "object { \n"
	    "  #include \"%s\" \n",objname);

    if (arg_texture->count > 0) {

      fprintf(scenefile,
	      "  texture {\n"
	      "     checker\n"
	      "     uv_mapping texture {\n"
	      "        pigment { Black }\n"
	      "        finish { Dull }\n"
	      "     }\n"
	      "     uv_mapping texture {\n"
	      "        pigment { White }\n"
	      "        finish { Dull }\n"
	      "     }\n"
	      "   }\n"
	      );
      
    } else { 

	fprintf(scenefile,
		"  texture {  \n"
		"     pigment { White } \n"
		"     finish { Dull } \n "
		"  }\n"
		"}\n");

      }

    fprintf(scenefile,"camera {\n");
    
    if (arg_print->count > 0) {
      fprintf(scenefile,"print_camera\n");
    } else if (arg_video->count > 0) {
      fprintf(scenefile,"video_camera\n");
    } else {
      fprintf(scenefile,"screen_camera\n");
    }

    fprintf(scenefile,"}\n\n");

    fprintf(scenefile,"background { rgb <1.0,1.0,1.0> }\n");
    fclose(scenefile);

    if (!QUIET) { printf("done.\n"); }

    /* We are now going to create a script to do the actual work. */

    if (!QUIET) { printf("povsnap: Building script to execute subcommands and POV-ray..."); }

    FILE *scriptfile;

    scriptfile = fopen_or_die("rerun","w");
    
    if (arg_texture->count > 0) {
      
      fprintf(scriptfile,
	      "cp ../%s .\n"
	      "cp ../%s .\n"
	      "/usr/local/bin/off2pov %s --uvfile\n"
	      "/usr/local/bin/makecam %s\n"
	      "/usr/local/bin/makelights %s\n",
	      basename,
	      mangle(basename,".off",".uv"),
	      basename,bboxname,bboxname);

    } else {
      
      if (!QUIET) {

       fprintf(scriptfile,
	      "cp ../%s .\n"
	      "/usr/local/bin/off2pov %s\n"
	      "/usr/local/bin/makecam %s\n"
	      "/usr/local/bin/makelights %s\n",
	      basename,
	      basename,bboxname,bboxname);


      } else {

       fprintf(scriptfile,
	      "cp ../%s .\n"
	      "/usr/local/bin/off2pov %s 1>/dev/null \n"
	      "/usr/local/bin/makecam %s 1>/dev/null \n"
	      "/usr/local/bin/makelights %s 1>/dev/null \n",
	      basename,
	      basename,bboxname,bboxname);

      }

    }

    fprintf(scriptfile,
	    "povray +FN %s ",scenename);

    if (arg_print->count > 0) {

      fprintf(scriptfile,"+W%d +H%d",
	      1024,768);

    } else if (arg_video->count > 0) {

      fprintf(scriptfile,"+W%d +H%d",
	      640,480);

    } else {

      fprintf(scriptfile,"+W%d +H%d",
	      320,240);

    }

    if (arg_display->count > 0) {

      fprintf(scriptfile," Display=Off ");

    }

    fprintf(scriptfile," +A All_file=povrayoutput.txt Verbose=true +o%s ",pngname);

    if (QUIET) {

      fprintf(scriptfile,"-GA 1>/dev/null 2>/dev/null \n");
      
    } else {

      fprintf(scriptfile,"+GA\n");

    }
    
    /* Eventually, we'll put some cropping stuff in here. */

    fclose(scriptfile);

    if (!QUIET) { printf("done.\n"); }

    /* Now execute the script... */

    if (!QUIET) { printf("povsnap: Calling programs and POV-ray. (Could take some time.) \n"); }

    chmod_or_die("rerun",S_IRWXU);
    system_or_die("./rerun");
    chdir_or_die("..");
    
    sprintf(cmdstring,"cp ./%s/%s .",dirname,pngname);
    system_or_die(cmdstring);
    
    if (QUIET) { printf("%s\n",pngname); } else {

      printf("File output to %s.\n",pngname);

    }

    /* If we generated a .tube.off file, clean it up here. */

  }

  exit(0);

}



