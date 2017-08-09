/*

	toros_ropelength: This procedure takes a (one-component) VECT
	file as input and computes the ropelength of the knot using TOROS.
	This code is only intended as a check for octrope... it won't run 
	unless testtoros is on the current path. 

*/

#ifdef HAVE_CONFIG_H
#include <../../config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_MATH_H
#include <math.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

#include<plCurve.h>
#include"../octropesrc/octrope.h"
#include<argtable2.h>

/* Global variables live here. */

struct arg_lit  *verbose;
struct arg_file *arg_infile;
struct arg_lit  *help;
struct arg_dbl  *arg_pass;

struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

FILE    *infile_fptr,*outfile_fptr;

int    QUIET=0;
double PASS = 1e-10;

double toros_thickness(plCurve *L)

{

  /* Step 0. Create a temporary file to save our guy in MING format. */

  char *template;

  template = calloc(32,sizeof(char));
  sprintf(template,"mingXXXXXX");
  int des;

  des = mkstemp(template);
  
  FILE *mingfile;
  mingfile = fdopen(des,"w");

  /* Step 1. Write the file in ming format. */

  fprintf(mingfile,"{%d}\n",L->cp[0].nv);

  int vt;
  for(vt=0;vt<L->cp[0].nv;vt++) {

    fprintf(mingfile,"  %21.17f  %21.17f  %21.17f\n",
	    plc_M_clist(L->cp[0].vt[vt]));

  }

  fclose(mingfile);

  /* Step 2. Call toros to deal with the problem. */

  FILE *Tpipe;
  char *toutname;

  toutname = tmpnam(NULL); /* Temporary file for TOROS output */
  char pipestr[1024];

  sprintf(pipestr,"./testtoros > %s",toutname);
  Tpipe = popen(pipestr,"w");

  assert(Tpipe != NULL);

  /* testtoros expects:

     filename
     dt
     q

  */

  fprintf(Tpipe,"%s\ndt\nq\n",template); 
  pclose(Tpipe);

  /* Step 3. Load and parse the results from TOROS. */

  FILE *toutput;

  toutput = fopen(toutname,"r");
  assert(toutput != NULL);

  double tthi = -1.0;
  char tline[4096];

  while (fgets(tline,4096,toutput)) {

    if (sscanf(tline,"Injectivity Radius: %lf",&tthi) == 1) {

      break;

    }

  }

  assert(tthi > 0);

  /* Step 4. Clean up all the stupid temporary files. */

  remove(toutname);
  remove(template);
  remove("boogie"); // This file is saved automatically by TOROS

  return tthi;
}


int main(int argc,char *argv[])
{
  int       infilenum,nerrors;
  plCurve   *link;
  bool      allpass = true;
  
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     arg_infile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     arg_pass = arg_dbl0("p","pass","<x>","pass test if TOROS and octrope differ by < x"),
     quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  fprintf(stderr,"torosropelength (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n\n");

  if (arg_nullcheck(argtable) != 0)
    printf("torosropelength: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"torosropelength");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("torosropelength checks the ropelength provided by TOROS against octrope.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  QUIET = quiet->count > 0; /* Register if we're in batch mode */

  if (arg_pass->count > 0) { PASS = arg_pass->dval[0]; }

  printf("Filename                         Verts  Octrope Thi    TOROS Thi    Diff < %5g \n",PASS);
  printf("-------------------------------------------------------------------------------\n");
  
  for(infilenum = 0;infilenum < arg_infile->count;infilenum++) {

    /* We start by loading the tube from core. */

    FILE *infile_fptr;

    infile_fptr = fopen(arg_infile->filename[infilenum],"r");
  
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"torosropelength: Couldn't open file %s.\n",
	      arg_infile->filename[infilenum]);
       continue;  /* Try the next file */
       
    }

    int plr_error_num;
    char plr_error_str[1024];
    
    link = plc_read(infile_fptr,
		    &plr_error_num,plr_error_str,sizeof(plr_error_str));
    
    /* We now demonstrate the octrope library's error handling protocol: */
  
    if (plr_error_num > 0) {   /* This is the signal for an error. */
      
      fprintf(stderr,"torosropelength: link reading error\n%s\n",plr_error_str);
      continue;  /* Try the next file */
      
    }
    
    fclose(infile_fptr);

    /* We now disqualify the file if it is open or contains > 1 component */

    if (link->nc > 1 || link->cp[0].open) { continue; }

    /* Report to user */

    printf("%-32s ",arg_infile->basename[infilenum]);

    double octThi;

    octThi = octrope_thickness(link,NULL,0,1.0);
    printf("%5d %12g",plc_num_verts(link),octThi); 

    /* Now we call TOROS to compare. */

    double torosThi;

    torosThi = toros_thickness(link);

    printf(" %12g",torosThi);

    if (fabs(octThi - torosThi) < PASS) { 

      printf("        pass\n");

    } else {

      printf("        FAIL\n");
      allpass = false;

    }

    plc_free(link);

  }

  printf("\nFinal Result: ");
  if (allpass) { printf("PASS\n"); exit(0); } else { printf("FAIL\n"); exit(1); }
  
}

