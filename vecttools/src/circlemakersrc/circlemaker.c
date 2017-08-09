/*******************************************************************

  circlemaker.c : Creates simple chains. 
	   
	   ************************************************************/

#include<circlemaker.h>
#include<../utilib/ordie.h>
#include<../utilib/mangle.h>
#include<argtable2.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];
extern int  VERBOSE;

void circlebuf(int verts,plc_vector center,double radius,plc_vector basis_a,plc_vector basis_b,plc_vector *buf)

/* Fills the buffer with a circle chosen carefully to have a reflection symmetry across basis_a. */

{
  double theta=0,t_step;
  double TWOPI = 6.2831853071795864769;

  t_step = TWOPI/(double)(verts);
  int i;

  for(i=0;i<verts;i++,theta+=t_step) {

    plc_M_vlincomb(buf[i],radius*cos(theta),basis_a,radius*sin(theta),basis_b);
    plc_M_add_vect(buf[i],center);

  }
}
  

plCurve *circle(int verts,double radius) 

{

  int cp; 
  
  int *nv;
  bool *open;
  int *cc;

  plCurve *tk;

  cp = 1;

  int i;

  nv = calloc_or_die(cp,sizeof(int));
  open = calloc_or_die(cp,sizeof(bool));
  cc = calloc_or_die(cp,sizeof(int));

  for (i=0;i<cp;i++) { 

    nv[i] = 2*(verts/2);  /* We need to make sure that verts is even to have reflection symmetry for odd chains */
    open[i] = false;
    cc[i] = 0;

  }

  tk = plc_new(cp,nv,open,cc);

  /* We are now prepared to build the actual link. The deal is this. */

  plc_vector basis_a[2] = {{{1,0,0}},{{0,1,0}}};
  
  circlebuf(nv[0],plc_build_vect(0,0,0),radius,basis_a[0],basis_a[1],tk->cp[0].vt);
  plc_fix_wrap(tk);

  free(nv);
  free(open);
  free(cc);

  return tk;

}

  
/**********************************************************************/

struct arg_lit  *verbose;
struct arg_file *outfile;
struct arg_lit  *help;

struct arg_dbl  *target_thickness;
struct arg_dbl  *radius;
struct arg_int  *verts;
struct arg_str  *axis;

struct arg_end  *end;
struct arg_end  *helpend;

plCurve *tk;
FILE    *outfile_fptr;
int     VERBOSE = 0;

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  char     outfile_name[1000];
  
  void *argtable[] = 
    {
      verbose = arg_lit0("v","verbose","print debugging information"),
      outfile = arg_file0("o","output","<filename>","output file name"),
      
      verts = arg_intn("v","verts","<int>",0,100,"verts"),  
      radius = arg_dbln("r","radius","<dbl>",0,100,"radius"),
      axis = arg_strn("a","axis","<x,y,z>",0,100,"normal axis"),
      
      help = arg_lit0("?","usage,help","display help message"),
      end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  int nerrors;
  
  plCurve *tknot;
  
  /* First, we parse the command-line arguments using argtable. */
  
  if (arg_nullcheck(argtable) != 0)
    printf("circlemaker: Insufficient memory to allocate argument table.\n");
  
  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"circlemaker");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("circlemaker generates collections of circles in coordinate planes\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (help->count > 0) {
    
    printf("circlemaker generates collections of circles in coordinate planes\n"
           "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);
    
  }
  
  double link_rad = 3;

  if (radius->count > 0) {
    
    link_rad = radius->dval[0];

  }

  int tkverts[100];
  int i;

  for(i=0;i<100;i++) { tkverts[i] = 150; }
  
  if (verts->count > 0) {

    for(i=0;i<verts->count;i++) {

      tkverts[i] = verts->ival[i];

    }

  }

  if (verbose->count > 0) {

    VERBOSE = 10;

  }

  printf("circlemaker\n"
	 "Generating circle with radius %g and %d vertices.\n",
	 link_rad,tkverts);

  tknot = circle(tkverts,link_rad); 

  /* Now we write the curve. */

  if (outfile->count > 0) {

    outfile_fptr = fopen(outfile->filename[0],"w");

    if (outfile_fptr == NULL) {

      fprintf(stderr,"circlemaker: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

    strcpy(outfile_name,outfile->filename[0]);


  } else {

    sprintf(outfile_name,"circle-%03g.vect",link_rad);
    outfile_fptr = fopen(outfile_name,"w");
    
    if (outfile_fptr == NULL) {
      
      fprintf(stderr,"circlemaker: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

  }
 
  plc_write(outfile_fptr,tknot);
  printf("Wrote circle to file %s\n\n",outfile_name);

  plc_free(tknot);
  fclose(outfile_fptr);

  exit(0);

}




