/*******************************************************************

 transformvect.c : Moves or rotates Geomview VECT files in space.
	   
	   ************************************************************/

#include<transformvect.h>
#include<argtable2.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];

void plc_translate(plCurve *curve,plc_vector v)

/* Translates the curve by v */

{
  int i,j;
  
  for(i=0;i<curve->nc;i++) {

    for(j=0;j<curve->cp[i].nv;j++) {

      curve->cp[i].vt[j] = plc_vect_sum(curve->cp[i].vt[j],v);

    }
    
  }

  plc_fix_wrap(curve);

}

/**********************************************************************/

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_str  *translate;
struct arg_end *end;
struct arg_end *helpend;

plCurve *core;
FILE         *infile_fptr,*outfile_fptr;

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,nerrors;
  char           outfile_name[1000],outfile_tail[1000];
  plc_vector     tvec;
 
  void *argtable[] = 
    {verbose = arg_lit0("v","verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     translate = arg_str0("t","translate","(x,y,z)","translates by vector (x,y,z)"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("transformtvect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"transformvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("transformvect moves or rotates a VECT curve in space\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }
  
  if (translate->count > 0) { 

    if (sscanf(translate->sval[0],"(%lf,%lf,%lf)",
	       &(tvec.c[0]),&(tvec.c[1]),&(tvec.c[2])) != 3) {

      printf("transformvect: Couldn't parse %s.\n",
	     translate->sval[0]);

      exit(1);

    }
 
  }
   
  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"transformvect: Couldn't open file %s.\n",
	       corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }
  
     int plr_error_num;
     char plr_error_str[1024];

     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"transformvect: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     /* Now we transform the link. */

     if (translate->count > 0) { 
       
       plc_translate(core,tvec);

     }
     
     /* Now we output */

     if (strlen(corefile->basename[infilenum]) > sizeof(outfile_name)-20) {
       
       fprintf(stderr,"transformvect: Ridiculously long input filename "
	       "can't be parsed.\n");
       exit(1);
       
     }
     
     sprintf(outfile_name,"%s",corefile->basename[infilenum]);
     
     if (strstr(outfile_name,".vect") != NULL) {
       
       sprintf(strstr(outfile_name,".vect"),".xform.vect");
       
     } else {
       
       sprintf(outfile_tail,".xform.vect");	 
       strcat(outfile_name,outfile_tail);
       
     }
     
     outfile_fptr = fopen(outfile_name,"w");
       
     if (outfile_fptr == NULL) {
       
       fprintf(stderr,
	       "transformvect: Couldn't open %s for writing.\n",outfile_name);
       exit(1);
       
     }
     
     /* Now write the core to the file. */
     
     plc_write(outfile_fptr,core);
     fclose(outfile_fptr);
     plc_free(core);
  
  }

  return 0;

}




