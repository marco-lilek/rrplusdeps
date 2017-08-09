
/*

  vect2maya : Converts a VECT file into a Maya MEL script which generates a NURBS curve.

*/

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

#ifdef HAVE_STDBOOL_H
#include<stdbool.h>
#endif

#include"plCurve.h"
#include<argtable2.h>

struct arg_lit  *verbose;
struct arg_file *corefile;
struct arg_file *outfile;
struct arg_lit  *help;
struct arg_dbl  *edgelength;
struct arg_dbl  *resolution;
struct arg_int  *verts;
struct arg_lit  *info;
struct arg_end  *end;
struct arg_end  *helpend;

plCurve *core;
FILE    *infile_fptr,*outfile_fptr;

int main(int argc,char *argv[])
{

  int   infilenum,nerrors;
   
  void *argtable[] = 
    {
     verbose = arg_lit0(NULL,"verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     outfile = arg_filen("o","output","<file>",0,100000,"output filename(s)"),
     // edgelength = arg_dbl0("e","edgelength","<x>","splines to edges of length < x"),
     // resolution = arg_dbl0("r","resolution","<n>","spline to n vertices per unit length"),
     // info = arg_lit0("i","info","display resolution information without splining"),
     // mps = arg_lit0(NULL,"mps","minrad preserving spline (instead of cubic)"),
     verts = arg_intn("v","verts","<n>",0,100000,"number of cvs (for each NURBS)"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  char outfile_name[16000],outfile_tok[16000], *tok;
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("vect2maya: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"vect2maya");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("vect2maya creates a MEL script which regenerates the curve as a NURBS\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);
    }
    
  }
  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"vect2maya: Couldn't open file %s.\n",corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }
  
     int plr_error_num;
     char plr_error_str[1024];

     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"vect2maya: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     printf("vect2maya: Read link from %s.\n",
	    corefile->basename[infilenum]);
   

     if (infilenum == 0) {

       sprintf(outfile_tok,"%s",corefile->basename[infilenum]);
       
       // We need to detect the .rN token if present in filename.
       // First, we delete any trailing .vect 
       
       tok = strstr(outfile_tok,".vect");
       if (tok != NULL) {*tok = 0;}  // Truncate the string, killing .vect
       
       // The string outfile_tok is now the truncated base of the new name
       
       sprintf(outfile_name,"%s.mel",outfile_tok);
       
       // We now scrap all this if the user has overridden the automatic filename generation
       
       if (outfile->count > 0) {
	 
	 sprintf(outfile_name,"%s",outfile->filename[infilenum%outfile->count]);
	 
       }
       
       outfile_fptr = fopen(outfile_name,"w");
       
       if (outfile_fptr == NULL) {
	 
	 fprintf(stderr,
		 "vect2maya: Couldn't open %s for writing.\n",outfile_name);
	 exit(1);
	 
       }

     }
       
     /* Now write each component to a file. */

     int cp;

     for(cp=0;cp < core->nc;cp++) {

       if (infilenum == 0) { /* This is the first curve in an animation set */

	 fprintf(outfile_fptr,"curve -degree 1 \n");
	 
	 int vt; 
	 for(vt=0;vt < core->cp[cp].nv; vt++) { 
	   
	   fprintf(outfile_fptr,"\t -p %g %g %g \n",plc_M_clist(core->cp[cp].vt[vt]));
	   
	 }
	 
	 fprintf(outfile_fptr," -n \"Component%d\" ;\n",cp);
	 
	 if (!core->cp[cp].open) { 

	   fprintf(outfile_fptr,"closeCurve -rpo \"Component%d\" ;\n",cp);

	 }

	 fprintf(outfile_fptr,"fitBspline -constructionHistory true -tolerance 0.01 -n \"SplineComponent%d\" \"Component%d\" ; \n",cp,cp);
	 fprintf(outfile_fptr,"\n");
	
       } else {

	 int vt; 
	 for(vt=0;vt < core->cp[cp].nv; vt++) { 
	   
	   fprintf(outfile_fptr,"setKeyframe -time %g -value << %g , %g , %g >> Component%d.cv[%d] ; \n",
		   infilenum/24.0,plc_M_clist(core->cp[cp].vt[vt]),cp,vt);

	 }
 
       }

     }
     
     plc_free(core);    
     
  }

  fclose(outfile_fptr);          

  if (infilenum == 1) {

    printf("vect2maya: Wrote link to %s.\n",
	   outfile_name);

  } else {

    printf("vect2maya: Animated collection of links to %s.\n",outfile_name);

  }

  if (verbose->count > 0) {

    printf("\nvect2maya processed %d input files.\n",infilenum);
    
  }

  return 0;

}

  
