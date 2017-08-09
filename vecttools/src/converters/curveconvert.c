/*******************************************************************

  curveconvert.c : Converts between file formats for polygonal knots.

	   ************************************************************/

#include<curveconvert.h>
#include<argtable2.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];

/**********************************************************************/

struct arg_lit  *verbose;
struct arg_file *corefile;
struct arg_file *outfile;
struct arg_str  *argfrom;
struct arg_str  *argto;
struct arg_lit  *help;
struct arg_end  *end;
struct arg_end  *helpend;

plCurve *core;
FILE    *infile_fptr,*outfile_fptr;
typedef enum  fileformat_type { SONO=0, MAYA=1, VECT=2, KNOTPLOT=3, KNT=4 } fileformat;
char    extensions[6][32] = { "txt", "maya", "vect", "coords", "knt" };

plCurve *plc_read_from_knotplot(FILE *infile)
{
  return NULL;
}

plCurve *plc_read_from_sono(FILE *infile)
{
  int nv;
  bool open = false;
  int cc = 0;
  int vt;
  plCurve *L;
  int c;

  if (fscanf(infile,"%d",&nv) != 1) {

    return NULL;

  } 

  L = plc_new(1,&nv,&open,&cc);

  for(vt=0;vt<nv;vt++) {

    for(c=0;c<3;c++) {

      if (fscanf(infile,"%lf",&(L->cp[0].vt[vt].c[c])) != 1) {

	plc_free(L);
	return NULL;

      }

    }

  }

  if (plc_distance(L->cp[0].vt[0],L->cp[0].vt[nv-1]) < 1e-4) {

    L->cp[0].nv--;  /* The first and last vertices are the same */

  }

  plc_fix_wrap(L);

  return L;

}

void plc_write_to_sono(FILE *outfile,plCurve *core)
{
  int vt,c;

  if (core->nc > 1) {

    fprintf(stderr,"curveconvert: WARNING! Sono output contains only the first component of %d component curve.\n",core->nc);

  }

  fprintf(outfile,"%d\n",core->cp[0].nv);
  
  for(vt=0;vt<core->cp[0].nv;vt++) {

    for(c=0;c<3;c++) {

      fprintf(outfile,"%23.17g\n",core->cp[0].vt[vt].c[c]);

    }

  }
    
}

void plc_write_to_maya(FILE *outfile,plCurve *core)
{
  int vt;

  if (core->nc > 1) {

    fprintf(stderr,"curveconvert: WARNING! Maya output contains only the first component of %d component curve.\n",core->nc);

  }
  
  for(vt=0;vt<core->cp[0].nv;vt++) {

    fprintf(outfile,"<< %23.17g, %23.17g, %23.17g >> \n",plc_M_clist(core->cp[0].vt[vt]));

  }

}

void plc_write_to_knotplot(FILE *outfile,plCurve *core)
{
}


/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,nerrors;
   
  void *argtable[] = 
    {
     verbose = arg_lit0(NULL,"verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     outfile = arg_filen("o","output","<file>",0,100000,"output filename(s)"),
     argfrom = arg_str1("f","from","<extension>","convert FROM format"),
     argto = arg_str0("t","to","<extension>","convert TO format"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("curveconvert: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"curveconvert");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("curveconvert converts VECT files to and from other formats\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  fileformat from,to;

  /* We now try to set the from and to variables from the command line (if possible). */

  if (argfrom->count > 0) {

    if (!strcmp(argfrom->sval[0],"VECT") || !strcmp(argfrom->sval[0],"vect")) {

      from = VECT;

    } else if (!strcmp(argfrom->sval[0],"SONO") || !strcmp(argfrom->sval[0],"sono")) {

      from = SONO;

    } else if (!strcmp(argfrom->sval[0],"KNT") || !strcmp(argfrom->sval[0],"knt")) {

      from = SONO;

    } else if (!strcmp(argfrom->sval[0],"KP") || !strcmp(argfrom->sval[0],"kp") || !strcmp(argfrom->sval[0],"KnotPlot") || !strcmp(argfrom->sval[0],"knotplot") || !strcmp(argfrom->sval[0],"coords") || !strcmp(argfrom->sval[0],"KNOTPLOT")) {

      from = KNOTPLOT;

    } else if (!strcmp(argfrom->sval[0],"MAYA") || !strcmp(argfrom->sval[0],"maya")) {

      from = MAYA;

    } else {

      fprintf(stderr,"curveconvert: Couldn't parse 'from' filetype %s.\n",argfrom->sval[0]);
      exit(1);

    }

  } else {

     fprintf(stderr,"curveconvert: Required to have a --from=<filetype>.\n");
     exit(1);

  }

  if (argto->count > 0) {

    if (!strcmp(argto->sval[0],"VECT") || !strcmp(argto->sval[0],"vect")) {

      to = VECT;

    } else if (!strcmp(argto->sval[0],"SONO") || !strcmp(argto->sval[0],"sono")) {

      to = SONO;

    } else if (!strcmp(argto->sval[0],"KP") || !strcmp(argto->sval[0],"kp") || !strcmp(argto->sval[0],"KnotPlot") || !strcmp(argto->sval[0],"knotplot") || !strcmp(argto->sval[0],"coords") || !strcmp(argto->sval[0],"KNOTPLOT")) {

      to = KNOTPLOT;

    } else if (!strcmp(argto->sval[0],"MAYA") || !strcmp(argto->sval[0],"maya")) {

      to = MAYA;

    } else {

      fprintf(stderr,"curveconvert: Couldn't parse 'to' filetype %s.\n",argto->sval[0]);
      exit(1);

    }

  } else {

    to = VECT;

  }

  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */
    
    infile_fptr = fopen(corefile->filename[infilenum],"r");
    
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"curveconvert: Couldn't open file %s.\n",corefile->filename[infilenum]);
      continue;  /* Try the next file */
      
    }
  
    int plr_error_num;
    char plr_error_str[1024];

    if (from == VECT) {
    
      core = plc_read(infile_fptr,
		      &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
      if (plr_error_num > 0) {   /* This is the signal for an error. */
	
	fprintf(stderr,"curveconvert: link reading error\n%s\n",plr_error_str);
	continue;  /* Try the next file */
	
      }

    } else if (from == KNOTPLOT) {

      core = plc_read_from_knotplot(infile_fptr);

    } else if (from == SONO) {

      core = plc_read_from_sono(infile_fptr);

    } else {

      core = NULL;

    }

    fclose(infile_fptr);

    if (core == NULL) {

      fprintf(stderr,"curveconvert: Could not load a file in format %s from %s.\n",extensions[from],corefile->basename[infilenum]);
      continue;

    }
  
    printf("curveconvert: Read link from %s.\n",
	   corefile->basename[infilenum]);
    
    /* Now we set the output name. */
    
    char outfile_name[16000],outfile_tok[16000], *tok;
    sprintf(outfile_tok,"%s",corefile->basename[infilenum]);
    
    // First, we delete any trailing extension corresponding to the old filename

    tok = strstr(outfile_tok,extensions[from]);
    if (tok != NULL) {*tok = 0;}  // Truncate the string.
    
    // The string outfile_tok is now the truncated base of the new name
    
    sprintf(outfile_name,"%s%s",outfile_tok,extensions[to]);
    
    // We now scrap all this if the user has overridden the automatic filename generation
    
    if (outfile->count > 0) {
      
      sprintf(outfile_name,"%s",outfile->filename[infilenum%outfile->count]);
      
    }
    
    outfile_fptr = fopen(outfile_name,"w");
    
    if (outfile_fptr == NULL) {
      
      fprintf(stderr,
	      "curveconvert: Couldn't open %s for writing.\n",outfile_name);
      exit(1);
      
    }
    
    /* We now do the actual conversion */
    
    if (to == SONO) {
      
      plc_write_to_sono(outfile_fptr,core);
      
    } else if (to == MAYA) {
      
      plc_write_to_maya(outfile_fptr,core);
      
    } else if (to == KNOTPLOT) {
      
      plc_write_to_knotplot(outfile_fptr,core);
      
    } else if (to == VECT) {
      
      plc_write(outfile_fptr,core);
      
    }
    
     fclose(outfile_fptr);
     
     printf("curveconvert: Wrote link to %s.\n",
	    outfile_name);
     plc_free(core);    
     
  }
  
  if (verbose->count > 0) {
    
    printf("\ncurveconvert processed %d input files.\n",infilenum);
    
  }
  
  return 0;
  
}




