/*******************************************************************

  claspmaker.c : Divides Geomview vect files into connected components,
                or splits them in other ways.
	   
	   ************************************************************/

#include<claspmaker.h>
#include<../utilib/ordie.h>
#include<../utilib/mangle.h>
#include<argtable2.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];
extern int  VERBOSE;

void ztranslate(plCurve *curve,int cp,double dist)

{ 
  int i;
  
  for(i=0;i<curve->cp[cp].nv;i++) {

    curve->cp[cp].vt[i].c[2] += dist;

  }

  plc_fix_wrap(curve);

}

plCurve *simpleclasp(double tail_length,
		     bool fixendpoints,bool planeendpoints) 

/* The simple clasp (with opening angle 90 degrees, or tau=1 and lambda=1) is a sufficiently */
/* special and weirdly degenerate case that we provide special code for it. */

{

  #define VERTS 100
  
  int arcnum;
  int nv[2] = {VERTS,VERTS};
  bool open[2] = {true,true};
  int cc[2] = {0,0};

  plCurve *clasp;

  clasp = plc_new(2,nv,open,cc);
  assert(clasp != NULL);

  plc_vector leftnormal[2],rightnormal[2];
  plc_vector buildx,buildy;
  double tailx,arcx;
  double taily,arcy;
  plc_vector left_endpoint,left_junction,
     right_endpoint,right_junction,arc_center;

  int tailverts,arcverts;
  double arclength,complength;

  double opening_angle = 1.4;
  double arcrad = 0.8;

  arclength = arcrad*2*opening_angle;
  complength = arclength + 2*tail_length;
  tailverts = (int)(floor(VERTS*tail_length/complength));
  arcverts = VERTS-2*tailverts;
 
  for(arcnum=0;arcnum<2;arcnum++) {
    
    if (arcnum == 0) {

      buildx = plc_build_vect(1,0,0);
      buildy = plc_build_vect(0,0,1);

    } else {

      buildx = plc_build_vect(0,1,0);
      buildy = plc_build_vect(0,0,-1);

    }
    
    tailx = 1.0;
    arcx  = sin(opening_angle)*arcrad; // Radius of the arc is now arcrad
    
    taily = sin(opening_angle)*tail_length;
    arcy = taily - cos(opening_angle)*arcrad;
    
    /* We have now located the most important points on the arc,
       the left endpoint (-tailx-arcx,0), the left junction (-arcx,taily),
       the right junction (arcx,taily), the right endpoint (tailx+arcx,0),
       and the center of the arc (0,arcy). */ 
    
    left_endpoint = plc_vlincomb(-tailx,buildx,0,buildy);
    left_junction = plc_vlincomb(-arcx,buildx,taily,buildy);
    right_junction = plc_vlincomb(arcx,buildx,taily,buildy);
    right_endpoint = plc_vlincomb(tailx,buildx,0,buildy);
    arc_center = plc_vlincomb(0,buildx,arcy,buildy);

    leftnormal[arcnum] =                      // The normals are always straight up and down.
      plc_normalize_vect(buildy,NULL);
    rightnormal[arcnum] = 
      plc_normalize_vect(buildy,NULL);
    
    /* Our job now is to fill in the vertices between these points. */
    
    int vt,i;
    
    for(i=0,vt=0;i<tailverts;i++,vt++) {
      
      clasp->cp[arcnum].vt[vt] = plc_vlincomb(1-(i/(double)(tailverts)),left_endpoint,
					      (i/(double)(tailverts)),left_junction);
    }
    
    double theta;
    
    for(i=0,theta=-opening_angle;i<arcverts;i++,theta+=2*opening_angle/(double)(arcverts),vt++) {
      
      clasp->cp[arcnum].vt[vt] = 
	plc_vect_sum(arc_center,
		     plc_vlincomb(arcrad*sin(theta),buildx,
				  arcrad*cos(theta),buildy));
      
    }
    
    for(i=1;i<tailverts+1;i++,vt++) {
      
      clasp->cp[arcnum].vt[vt] = plc_vlincomb(1-(i/(double)(tailverts)),right_junction,
					      i/(double)(tailverts),right_endpoint);
      
    }

  }

  plc_fix_wrap(clasp);

  /* For this very special configuration, we just compute where we 
     want the second component to go. */

  ztranslate(clasp,1,2*arcy+arcrad-0.05);

  if (VERBOSE > 5) {

    printf("Now adjusting z position of second component...\n");
    printf("Thickness of final configuration is %g.\n",octrope_thickness(clasp,NULL,0,1));
  }
  
  /* we now need to set the constraints appropriately */

  if (fixendpoints) {

    plc_set_fixed(clasp,0,0,clasp->cp[0].vt[0]);
    plc_set_fixed(clasp,0,clasp->cp[0].nv-1,clasp->cp[0].vt[clasp->cp[0].nv-1]);

    plc_set_fixed(clasp,1,0,clasp->cp[1].vt[0]);
    plc_set_fixed(clasp,1,clasp->cp[0].nv-1,clasp->cp[1].vt[clasp->cp[0].nv-1]);

  } 

  if (planeendpoints) {

    plc_constrain_to_plane(clasp,0,0,1,leftnormal[0],
			   plc_dot_prod(leftnormal[0],clasp->cp[0].vt[0]));
    plc_constrain_to_plane(clasp,0,clasp->cp[0].nv-1,1,rightnormal[0],
			   plc_dot_prod(rightnormal[0],clasp->cp[0].vt[clasp->cp[0].nv-1]));
    
    plc_constrain_to_plane(clasp,1,0,1,leftnormal[1],
			   plc_dot_prod(leftnormal[1],clasp->cp[1].vt[0]));
    plc_constrain_to_plane(clasp,1,clasp->cp[1].nv-1,1,rightnormal[1],
			   plc_dot_prod(rightnormal[1],clasp->cp[1].vt[clasp->cp[1].nv]));
  
  }
    
    
  return clasp;

}

    
plCurve *clasp(double tail_length,double opening_angle,
	       double twist_angle,double lambda,
	       bool fixendpoints,bool planeendpoints) 

/* Build a clasp with given twist angle, tail length, and opening
   angle, with arc having radius lambda. */

{

  int arcnum;
  int nv[2] = {300,300};
  bool open[2] = {true,true};
  int cc[2] = {0,0};

  plCurve *clasp;

  clasp = plc_new(2,nv,open,cc);
  assert(clasp != NULL);

  plc_vector leftnormal[2],rightnormal[2];
  plc_vector buildx,buildy;
  double tailx,arcx;
  double taily,arcy;
  plc_vector left_endpoint,left_junction,
     right_endpoint,right_junction,arc_center;

  for(arcnum=0;arcnum<2;arcnum++) {
    
    if (arcnum == 0) {

      buildx = plc_build_vect(1,0,0);
      buildy = plc_build_vect(0,0,1);

    } else {

      buildx = plc_build_vect(cos(twist_angle),sin(twist_angle),0);
      buildy = plc_build_vect(0,0,-1);

    }
    
    tailx = cos(opening_angle)*tail_length;
    arcx  = sin(opening_angle)*lambda;
    
    taily = sin(opening_angle)*tail_length;
    arcy = taily - cos(opening_angle)*lambda;
    
    /* We have now located the most important points on the arc,
       the left endpoint (-tailx-arcx,0), the left junction (-arcx,taily),
       the right junction (arcx,taily), the right endpoint (tailx+arcx,0),
       and the center of the arc (0,arcy). */ 
    
    left_endpoint = plc_vlincomb(-tailx-arcx,buildx,0,buildy);
    left_junction = plc_vlincomb(-arcx,buildx,taily,buildy);
    right_junction = plc_vlincomb(arcx,buildx,taily,buildy);
    right_endpoint = plc_vlincomb(tailx+arcx,buildx,0,buildy);
    arc_center = plc_vlincomb(0,buildx,arcy,buildy);

    leftnormal[arcnum] = 
      plc_normalize_vect(plc_vect_diff(left_junction,left_endpoint),NULL);
    rightnormal[arcnum] = 
      plc_normalize_vect(plc_vect_diff(right_junction,right_endpoint),NULL);
    
    /* Our job now is to fill in the vertices between these points. */
    
    int vt,i;
    
    for(i=0,vt=0;i<100;i++,vt++) {
      
      clasp->cp[arcnum].vt[vt] = plc_vlincomb(1-(i/100.0),left_endpoint,
					      (i/100.0),left_junction);
    }
    
    double theta;
    
    for(i=0,theta=-opening_angle;i<100;i++,theta+=2*opening_angle/100.0,vt++) {
      
      clasp->cp[arcnum].vt[vt] = 
	plc_vect_sum(arc_center,
		     plc_vlincomb(lambda*sin(theta),buildx,
				  lambda*cos(theta),buildy));
      
    }
    
    for(i=1;i<101;i++,vt++) {
      
      clasp->cp[arcnum].vt[vt] = plc_vlincomb(1-(i/100.0),right_junction,
					      i/100.0,right_endpoint);
      
    }

  }

  plc_fix_wrap(clasp);

  /* Now we need to move the second component up in order to get the
     arcs to touch correctly in the center. It's not completely
     obvious how far to move up the curve, so we simply edge it up
     until the thickness starts to drop */

  if (VERBOSE > 5) {

    printf("Now adjusting z position of second component...\n");

  }

  double zmin=0,zmax=2*arcy + 2*lambda-0.5;

  for(;zmax - zmin > 0.01;) {

    double a,b,c;
    double zmid;

    zmid = (zmax+zmin)/2.0;

    ztranslate(clasp,1,zmin);
    a = octrope_thickness(clasp,NULL,0,lambda);

    ztranslate(clasp,1,(zmax-zmin)/2.0);
    b = octrope_thickness(clasp,NULL,0,lambda);

    ztranslate(clasp,1,(zmax-zmin)/2.0);
    c = octrope_thickness(clasp,NULL,0,lambda);

    ztranslate(clasp,1,-zmax);
    
    if (b < 0.51) { zmax = zmid; } else { zmin = zmid;}
  
  }
  
  ztranslate(clasp,1,zmin);
  
  /* we now need to set the constraints appropriately */

  if (fixendpoints) {

    plc_set_fixed(clasp,0,0,clasp->cp[0].vt[0]);
    plc_set_fixed(clasp,0,299,clasp->cp[0].vt[299]);

    plc_set_fixed(clasp,1,0,clasp->cp[1].vt[0]);
    plc_set_fixed(clasp,1,299,clasp->cp[1].vt[299]);

  } 

  if (planeendpoints) {

    plc_constrain_to_plane(clasp,0,0,1,leftnormal[0],
			   plc_dot_prod(leftnormal[0],clasp->cp[0].vt[0]));
    plc_constrain_to_plane(clasp,0,299,1,rightnormal[0],
			   plc_dot_prod(rightnormal[0],clasp->cp[0].vt[299]));
    
    plc_constrain_to_plane(clasp,1,0,1,leftnormal[1],
			   plc_dot_prod(leftnormal[1],clasp->cp[1].vt[0]));
    plc_constrain_to_plane(clasp,1,299,1,rightnormal[1],
			   plc_dot_prod(rightnormal[1],clasp->cp[1].vt[299]));
  
  }
    
    
  return clasp;

}

/**********************************************************************/

struct arg_lit  *verbose;
struct arg_file *outfile;
struct arg_lit  *help;

struct arg_dbl  *tau;
struct arg_dbl  *opening_angle;
struct arg_dbl  *twist_angle;
struct arg_dbl  *lambda;
struct arg_dbl  *tail_len;

struct arg_end  *end;
struct arg_end  *helpend;

struct arg_lit  *fixed_endpoints;
struct arg_lit  *planar_endpoints;
struct arg_lit  *simple_clasp;

plCurve *claspCurve;
FILE    *outfile_fptr;
int     VERBOSE = 0;

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  double   PI = 3.14159265358973;
  char     outfile_name[1000];
  double   openingAngle = {PI/2.0},twistAngle={PI/2.0};
  double   Lambda={1.0},tailLength = {2.0};
  bool     fixedEndpoints = {true}, planarEndpoints = {false};
 
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     outfile = arg_file0("o","output","<filename>","output file name"),
     simple_clasp = arg_lit0("s","simpleclasp","construct \"simple\" clasp with tau=1, lambda=1"),
     tau = arg_dbl0("t","tau","<0..1>","sin of opening angle for clasp"),
     opening_angle = arg_dbl0(NULL,"openingangle","<0..2pi>","opening angle for clasp"),
     twist_angle = arg_dbl0(NULL,"twistangle","<0..2pi>","twist components relative to one another by this angle"),
     lambda = arg_dbl0("l","lambda",">=1.0","stiffness for clasp"),
     tail_len = arg_dbl0(NULL,"taillength",">=0","length of straight tails for clasp"),
     fixed_endpoints = arg_lit0("f","fixed","endpoints fixed"),
     planar_endpoints = arg_lit0("p","planar","endpoints constrained to planes"),
     help = arg_lit0("?","usage,help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  int nerrors;

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("claspmaker: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"claspmaker");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("claspmaker generates simple clasp curves with "
	     "specified geometry\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (help->count > 0) {

    printf("claspmaker generates simple clasp curves with "
	   "specified geometry\n"
           "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);
    
  }

  if (tau->count > 0) {

    openingAngle = asin(tau->dval[0]);

  }

  if (opening_angle->count > 0) {

    openingAngle = opening_angle->dval[0];

  }

  if (twist_angle->count > 0) {

    twistAngle = twist_angle->dval[0];

  }

  if (lambda->count > 0) {

    Lambda = lambda->dval[0];

  }

  if (tail_len->count > 0) {

    tailLength = tail_len->dval[0];

  }

  if (fixed_endpoints->count > 0) {

    fixedEndpoints = true;
    planarEndpoints = false;

  }

  if (planar_endpoints->count > 0) {

    planarEndpoints = true;
    fixedEndpoints = false;

  }

  if (verbose->count > 0) {

    VERBOSE = 10;

  }

  if (simple_clasp->count > 0) {

    printf("claspmaker\n"
	   "Generating simple clasp with \n");

  } else {

    printf("claspmaker\n"
	   "Generating clasp with \n");

  }

  printf("\t tau    = %g\n"
	 "\t lambda = %g\n"
	 "\t twist  = %g\n"
	 "\t tail   = %g\n"
	 "\t endpts = ",
	 sin(openingAngle),Lambda,twistAngle,tailLength);
  
  if (fixedEndpoints) { printf("fixed.\n\n"); } 
  else if (planarEndpoints) { printf("planar.\n\n"); }
  else {printf("unknown. Error!\n\n"); exit(1); }
  
  if (simple_clasp->count > 0) {

    claspCurve = simpleclasp(tailLength,fixedEndpoints,planarEndpoints);

  } else {

    claspCurve = clasp(tailLength,openingAngle,twistAngle,
		       Lambda,fixedEndpoints,planarEndpoints);

  }
    
  if (outfile->count > 0) {

    outfile_fptr = fopen(outfile->filename[0],"w");

    if (outfile_fptr == NULL) {

      fprintf(stderr,"claspmaker: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

    strcpy(outfile_name,outfile->filename[0]);


  } else {

    sprintf(outfile_name,"clasp-tau-%4.3f-lambda-%4.3f-twist-%4.3f.vect",
	    sin(openingAngle),Lambda,twistAngle);
    
    outfile_fptr = fopen(outfile_name,"w");
    
    if (outfile_fptr == NULL) {
      
      fprintf(stderr,"claspmaker: could not open %s for writing.\n",
	      outfile->filename[0]);
      exit(1);

    }

  }
 
  plc_write(outfile_fptr,claspCurve);
  printf("Wrote clasp to file %s\n\n",outfile_name);


  plc_free(claspCurve);
  fclose(outfile_fptr);

  exit(0);

}




