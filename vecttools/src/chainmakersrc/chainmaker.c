/*******************************************************************

  chainmaker.c : Creates simple chains. 
	   
	   ************************************************************/

#include<chainmaker.h>
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
  

plCurve *chain(int links,int verts,bool waist,double link_radius) 

{

  int cp; 
  
  int *nv;
  bool *open;
  int *cc;

  plCurve *tk;

  cp = ((waist) ? 1 : 0);
  cp += links;

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

  plc_vector basis_b = {{0,0,1}}, basis_a[2] = {{{1,0,0}},{{0,1,0}}};
  double center_z;
  int bat = 0;
  
  for(i=0,center_z = -(1.5)*(link_radius)*(links-1)/2.0;i<links;i++,center_z += (1.5)*link_radius,bat = 1 - bat) {
  
    circlebuf(nv[i],plc_build_vect(0,0,center_z),link_radius,basis_a[bat],basis_b,tk->cp[i].vt);

  }

  /* Now if we're to put in a waist circle, it is always in the x-y plane. */

  if (waist) { 

    circlebuf(nv[i],plc_build_vect(0,0,0),5,basis_a[0],basis_a[1],tk->cp[i].vt);

  }

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
struct arg_dbl  *link_radius;
struct arg_int  *links;
struct arg_lit  *waist;
struct arg_int  *verts;

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
      
      links = arg_int1("l","links","<int>","number of links in the chain"),
      waist = arg_lit0("w","waist","add a horizontal 'waist' circle around the chain"),

      verts = arg_int0("v","verts","<int>","verts per component"),
      
      link_radius = arg_dbl0("r","linkradius","<dbl>","radius of links"),
      target_thickness = arg_dbl0("t","thickness","<dbl>","target thickness"),
      
      help = arg_lit0("?","usage,help","display help message"),
      end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  int nerrors;
  
  plCurve *tknot;
  
  /* First, we parse the command-line arguments using argtable. */
  
  if (arg_nullcheck(argtable) != 0)
    printf("chainmaker: Insufficient memory to allocate argument table.\n");
  
  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"chainmaker");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("chainmaker generates simple chains of a variable number of links with specified geometry\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (help->count > 0) {
    
    printf("chainmaker generates simple chains of a variable number of links with specified geometry\n"
           "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);
    
  }
  
  double target_thi = 0.75;
  
  if (target_thickness->count > 0) { 
    
    target_thi = target_thickness->dval[0];
    
  }
  
  double link_rad = 3;

  if (link_radius->count > 0) {
    
    link_rad = link_radius->dval[0];

  }

  int tkverts = 150;

  if (verts->count > 0) {

    tkverts = verts->ival[0];

  }

  if (verbose->count > 0) {

    VERBOSE = 10;

  }

  if (links->ival[0] < 1) {

    printf("chainmaker: links must be >= 1\n");

  }

  printf("chainmaker\n"
	 "Generating simple chain with \n");
  
  printf("\t links     = %d\n"
	 "\t linkradius= %g\n"
	 "\t verts     = %d\n"
	 "\t thickness = %g\n",
	 links->ival[0],link_rad,tkverts,target_thi);

  if (waist->count > 0) {

    printf("\t waist     = yes\n");

  } else {

    printf("\t waist     = no\n");

  }

  tknot = chain(links->ival[0],tkverts,(waist->count > 0),link_rad); 

  double thi;
  thi = octrope_thickness(tknot,NULL,0,1.0);

  plc_scale(tknot,target_thi/thi);

  /* Now we write the curve. */

  if (outfile->count > 0) {

    outfile_fptr = fopen(outfile->filename[0],"w");

    if (outfile_fptr == NULL) {

      fprintf(stderr,"chainmaker: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

    strcpy(outfile_name,outfile->filename[0]);


  } else {

    if (waist->count > 0) {
      
      sprintf(outfile_name,"%d-link-chain-waist.vect",
	      links->ival[0]);

    } else {

      sprintf(outfile_name,"%d-link-chain.vect",
	      links->ival[0]);

    }
    
    outfile_fptr = fopen(outfile_name,"w");
    
    if (outfile_fptr == NULL) {
      
      fprintf(stderr,"chainmaker: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

  }
 
  plc_write(outfile_fptr,tknot);
  printf("Wrote simple chain to file %s\n\n",outfile_name);

  plc_free(tknot);
  fclose(outfile_fptr);

  exit(0);

}




