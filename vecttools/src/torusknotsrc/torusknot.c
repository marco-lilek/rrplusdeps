/*******************************************************************

  torusknot.c : Creates (p,q) torus knots. 
	   
	   ************************************************************/

#include<torusknot.h>
#include<../utilib/ordie.h>
#include<../utilib/mangle.h>
#include<argtable2.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];
extern int  VERBOSE;

int gcd(int x, int y)
{
  int t;

  while (y) {
    t = x;
    x = y;
    y = t % y;
  }
  return(x);
}

plCurve *torusknot(int verts,int p,int q,double major_radius,double minor_radius,int symmetry) 

{

  int cp; 
  int vt;
  
  int *nv;
  bool *open;
  int *cc;

  plCurve *tk;

  cp = gcd(p,q);
  
  int pcheck;
  pcheck = gcd(symmetry,p);
  if (pcheck != symmetry && symmetry > 0) { 
    fprintf(stderr,"torusknot: Error. Desired symmetry %d must divide p = %d.\n",symmetry,p);
    exit(1);
  }

  int i;

  nv = calloc_or_die(cp,sizeof(int));
  open = calloc_or_die(cp,sizeof(bool));
  cc = calloc_or_die(cp,sizeof(int));

  for (i=0;i<cp;i++) { 

    nv[i] = p*q*(verts/(p*q));  /* We need to make sure that verts is a multiple of pq to have symmetry. */
    open[i] = false;
    cc[i] = 0;

  }

  tk = plc_new(cp,nv,open,cc);

  /* We are now prepared to build the actual knot */

  double pofs;
  double theta; 
  double tstep;
  double pi = 3.1415926535897932384626433;
  double pangle,qangle;
  double symmetryrad = 1;

  for (i=0;i<cp;i++) {

    for(vt=0,theta=0,pofs=i*(2*pi/q),tstep = 2*pi/(cp*(double)(nv[i]));
	vt<nv[i];
	vt++,theta+=tstep) {

      plc_vector loc;

      pangle = pofs + p*theta;
      qangle = q*theta;

      if (symmetry == q || symmetry == 0) { /* Leave things alone */
	symmetryrad = 1.0;
      } else {
	symmetryrad = 1.0 + 0.2*sin(symmetry*qangle);
      }
	
      loc = plc_build_vect(
			   symmetryrad*major_radius*cos(qangle)*(1+(minor_radius/major_radius)*cos(pangle)),
			   symmetryrad*major_radius*sin(qangle)*(1+(minor_radius/major_radius)*cos(pangle)),
			   minor_radius*sin(pangle)
			   );

      tk->cp[i].vt[vt] = loc;

    }
    
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
struct arg_dbl  *aspect_ratio;
struct arg_int  *q;
struct arg_int  *p;
struct arg_int  *s;
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
      
      p = arg_int1("p","pvalue","<int>","number of times around the minor radius"),
      q = arg_int1("q","qvalue","<int>","number of times around the major radius"),
      s = arg_int0("s","symmetry","<int>","rotational symmetry around major radius (must divide q)"),

      verts = arg_int0("v","verts","<int>","verts per component"),
      
      aspect_ratio = arg_dbl0("a","aspectratio","<dbl>","ratio of major/minor radii"),
      target_thickness = arg_dbl0("t","thickness","<dbl>","target thickness"),
      
      help = arg_lit0("?","usage,help","display help message"),
      end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  int nerrors;
  
  plCurve *tknot;
  
  /* First, we parse the command-line arguments using argtable. */
  
  if (arg_nullcheck(argtable) != 0)
    printf("torusknot: Insufficient memory to allocate argument table.\n");
  
  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"torusknot");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("torusknot generates torus knots with specified geometry\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (help->count > 0) {
    
    printf("torusknot generates torus knots with specified geometry and symmetry\n"
           "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);
    
  }

  int symmetry=0;

  if (s->count > 0) {

    symmetry = s->ival[0];

  }
  
  double target_thi = 0.75;
  
  if (target_thickness->count > 0) { 
    
    target_thi = target_thickness->dval[0];
    
  }
  
  double aspect_rat = 3;

  if (aspect_ratio->count > 0) {
    
    aspect_rat = aspect_ratio->dval[0];

  }

  int tkverts = 200;

  if (verts->count > 0) {

    tkverts = verts->ival[0];

  }

  if (verbose->count > 0) {

    VERBOSE = 10;

  }

  if (p->ival[0] < 1 || q->ival[0] < 1) {

    printf("torusknot: p and q must be >= 1\n");

  }

  if (tkverts < 6*p->ival[0]) { tkverts = 6*p->ival[0]; }
  if (tkverts < 6*q->ival[0]) { tkverts = 6*q->ival[0]; }

  printf("torusknot\n"
	 "Generating torus knot with \n");

  
  printf("\t p         = %d\n"
	 "\t q         = %d\n"
	 "\t aspect    = %g\n"
	 "\t verts     = %d\n"
	 "\t thickness = %g\n"
	 "\t symmetry  = %d-fold\n",
	 p->ival[0],q->ival[0],aspect_rat,tkverts,target_thi,symmetry);
  
  tknot = torusknot(tkverts,p->ival[0],q->ival[0],aspect_rat,1.0,symmetry);

  double thi;
  thi = octrope_thickness(tknot,NULL,0,1.0);

  plc_scale(tknot,target_thi/thi);

  printf("\t major    = %g\n"
	 "\t minor    = %g\n",
	 aspect_rat*target_thi/thi,target_thi/thi);

  /* Now we write the curve. */

  if (outfile->count > 0) {

    outfile_fptr = fopen(outfile->filename[0],"w");

    if (outfile_fptr == NULL) {

      fprintf(stderr,"torusknot: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

    strcpy(outfile_name,outfile->filename[0]);


  } else {

    sprintf(outfile_name,"torusknot-%d-%d.vect",
	    p->ival[0],q->ival[0]);
    
    outfile_fptr = fopen(outfile_name,"w");
    
    if (outfile_fptr == NULL) {
      
      fprintf(stderr,"torusknot: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

  }
 
  plc_write(outfile_fptr,tknot);
  printf("Wrote torus knot to file %s\n\n",outfile_name);

  plc_free(tknot);
  fclose(outfile_fptr);

  exit(0);

}




