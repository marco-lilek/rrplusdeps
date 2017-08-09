/*******************************************************************

perturbvect.c : Randomly varies vertices of plCurve. 
	   
	   ************************************************************/

#include<perturbvect.h>
#include<argtable2.h>
#include<../utilib/mangle.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];

/**********************************************************************/

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_dbl  *arg_radius;
struct arg_lit  *arg_normal;
struct arg_end *end;
struct arg_end *helpend;

plCurve  *core;
FILE     *infile_fptr,*outfile_fptr;
double    radius = 0.01;
bool      normal = false;

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,nerrors;
   
  void *argtable[] = 
    {verbose = arg_lit0(NULL,"verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     arg_radius = arg_dbl0("r","radius","<x>","radius of perturbation"),
     arg_normal = arg_lit0("n","normal","perturb in normal disks to each point"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("perturbvect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"perturbvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("perturbvect randomly varies the unconstrained vertices of a plCurve\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (arg_radius->count > 0) { radius = arg_radius->dval[0]; }
  if (arg_normal->count > 0) { normal = true; }
  
  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"constrainvect: Couldn't open file %s.\n",
	       corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }
  
     int plr_error_num;
     char plr_error_str[1024];

     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"constrainvect: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     /* Now we work on the link. */

     plc_perturb(core,radius);

     /* We have applied all constraints. Save the file. */

     printf("perturbvect: perturbed %s by radius of %g\n",
	    corefile->basename[infilenum],radius);

     char *newname;
     newname = mangle(corefile->filename[infilenum],
		      ".vect",".perturbed.vect");
     
     outfile_fptr = fmangle(corefile->filename[infilenum],
			    ".vect",".perturbed.vect");
     
     plc_write(outfile_fptr,core);

     printf("constrainvect: Saved file to %s.\n",newname);

     free(newname);
     plc_free(core);
     
  }

  return 0;

}




