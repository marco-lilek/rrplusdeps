/*******************************************************************

  smoothvect.c : Smooths vect files, being careful to reduce ropelength.

	   ************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#include<plCurve.h>
#include<octrope.h>
#include<../utilib/mangle.h>
#include<../utilib/ordie.h>
#include<argtable2.h>

/* Global variables live here. */

struct arg_int *steps;
struct arg_dbl *lambda;

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_file *outfile;
struct arg_lit  *help;
struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

plCurve *core;
FILE    *infile_fptr,*outfile_fptr;

int    QUIET=0;
double gLambda = 1.0;
int    maxSteps = 1000;


/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int         infilenum,nerrors;
  
  void *argtable[] = 
    {
      steps = arg_int0("s","steps","<int>","number of steps of smoothing to run"),
      lambda = arg_dbl0("l","lambda","<dbl>","stiffness parameter (1.0 for standard rope)"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
      outfile = arg_filen("o","outfile","<file>",0,100000,"output filenames"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("smoothvect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"smoothvect (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"smoothvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"smoothvect (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      printf("smoothvect smooths corners in vect file while decreasing ropelength.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  QUIET = quiet->count > 0; /* Register if we're in batch mode */
  if (!QUIET) { fprintf(stderr,"smoothvect (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n"); }

  if (lambda->count > 0) {

    gLambda = lambda->dval[0];

  }
 
  if (steps->count > 0) { maxSteps = steps->ival[0]; }

  /* Now we have parsed the arguments and are ready to work. */

  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"smoothvect: Couldn't open file %s.\n",
	       corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }

     int plr_error_num;
     char plr_error_str[1024];
  
     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     /* We now demonstrate the plCurve library's error handling protocol: */
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"smoothvect: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     /* Report to the user */

     if (!QUIET) {

       printf("%s\n",corefile->basename[infilenum]);

     }

     /* Now we actually do the smoothing. */

     double rop,thi,len,minrad,minstrut,lastrop;
     octrope_mrloc min_rad_locs[10];
     int n_minradlocs;
     int stepcount;
       
     octrope(core,&rop,&thi,&len,&minrad,&minstrut,
	     0,0,min_rad_locs,10,&n_minradlocs,
	     0,0,NULL,0,NULL,
	     NULL,0,gLambda);
     
     printf("Starting smooth with minrad = %g, minstrut = %g, rop = %g.\n",minrad,minstrut,rop);
       
     for(lastrop=rop,stepcount = 0;stepcount < maxSteps && lastrop >= rop;stepcount++) {

       lastrop = rop;

       int rollfwd;

       for(rollfwd = 0;rollfwd < 50;rollfwd++) {

	 octrope(core,&rop,&thi,&len,&minrad,&minstrut,
		 0,0,min_rad_locs,10,&n_minradlocs,
		 0,0,NULL,0,NULL,
		 NULL,0,gLambda);
	 
	 printf("%03.4g %03.4g %03.4g\n",rop,minrad,minstrut);
	 
	 //if (gLambda*minrad > 0.3*minstrut) { /* minrad is in control, but not by much */ break; }
	 
	 /* We are now actually going to smooth. We go to the vertex in control of minrad,
	    and blend it with the midpoint of the adjacent vertices. */
	 
	 core->cp[min_rad_locs[0].component].vt[min_rad_locs[0].vert] = 
	   plc_vweighted(0.8,
			 core->cp[min_rad_locs[0].component].vt[min_rad_locs[0].vert],
			 plc_vweighted(0.5,
				       core->cp[min_rad_locs[0].component].vt[min_rad_locs[0].vert+1],
				       core->cp[min_rad_locs[0].component].vt[min_rad_locs[0].vert-1])
			 );
	 
	 plc_fix_wrap(core);  /* Having changed a vertex, we fixwrap. */

       }
	 
     }


     if (!QUIET) {

       printf("%d steps of smoothing...\n",stepcount);

       octrope(core,&rop,&thi,&len,&minrad,&minstrut,
	       0,0,min_rad_locs,10,&n_minradlocs,
	       0,0,NULL,0,NULL,
	       NULL,0,gLambda);

       printf("Ending smooth with minrad = %g, minstrut = %g, rop = %g.\n",minrad,minstrut,rop);
       
     }	     
	 
     /* Now we output the smoothed version. */

     const char *outfilename;

     if (outfile->count > 0) {

       outfilename = outfile->basename[infilenum % outfile->count];

     } else {

       outfilename = mangle(corefile->basename[infilenum],"vect","smooth.vect");

     }

     FILE *outfile;
     outfile = fopen_or_die(outfilename,"w");
     plc_write(outfile,core);
     fclose(outfile);

     if (!QUIET) {

       printf("Wrote smoothed file to %s.\n",outfilename);

     }

     plc_free(core);

  }
  
  return 0;

}




