/*******************************************************************

  tube.c : Draws tubular neighborhoods of Geomview VECT curves.
           See the man page for details.
	   

	   ************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<tube.h>
#include<uv.h>
#include<../utilib/mangle.h>
#include<../utilib/ordie.h>
#include<argtable2.h>
#include<maxmin.h>

/* Global variables live here. */

struct arg_dbl *r;
struct arg_dbl *g;
struct arg_dbl *b;

struct arg_dbl *ratio;
struct arg_int *minsides;

struct arg_lit *capped;
struct arg_lit *closed;
struct arg_lit *equalized;
struct arg_lit *verbose;
struct arg_lit *makerings;

struct arg_lit *flared;
struct arg_dbl *flare_power;
struct arg_dbl *flare_radius;
struct arg_dbl *flare_distance;

struct arg_lit *bulged;
struct arg_dbl *bulge_power;
struct arg_dbl *bulge_radius;
struct arg_dbl *bulge_start;
struct arg_dbl *bulge_end;

struct arg_dbl  *radius;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_lit  *uv;

struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

plCurve *core;
FILE         *infile_fptr,*outfile_fptr;

int    nbulges;
double *bstarts,*bends,*bpowers,*brads;

double DEFAULT_BSTART = 1.0/3.0;
double DEFAULT_BEND   = 2.0/3.0;
double DEFAULT_BPOWER = 3.0;
double DEFAULT_BRAD   = 2.0;

int    nflares;
double *fpowers, *frads, *fdists;

double DEFAULT_FPOWER = 3.0;
double DEFAULT_FRAD   = 2.0;
double DEFAULT_FDIST  = 0.1;

int    QUIET=0;
int    MINSIDES=6;

int    MAKERINGS=0;
plCurve *rings = NULL;

void parse_color_args(plCurve L)

     /* This procedure parses the various occurences of the color 
	arguments and assigns colors to the components of L. We use
	the rule that the first occurence of a color refers to the 
	first component of the link and so forth... */

     /* If a color is specified, it overrides whatever color was 
	present in the input link. Otherwise, we leave that color
	information in place. */

{ 
  int ncolors;
  int i;

  ncolors = intmax(3,r->count,g->count,b->count);

  for(i=0;i<L.nc && i<ncolors;i++) {

    /* Make sure there is space for the new color info */

    if (L.cp[i].cc == 0) {

      L.cp[i].clr = (plc_color *)(malloc(1*sizeof(plc_color)));
      assert(L.cp[i].clr != NULL);

    }

    L.cp[i].cc = 1;

    /* Now assign R, G, and B values (if present). */

    L.cp[i].clr[0].r = (i < r->count) ? (r->dval[i]) : 1.0;
    L.cp[i].clr[0].g = (i < g->count) ? (g->dval[i]) : 1.0;
    L.cp[i].clr[0].b = (i < b->count) ? (b->dval[i]) : 1.0;
    L.cp[i].clr[0].alpha = 1.0;

  }

}

void parse_bulge_and_flare_args()

/* This procedure parses the bulge and flare arguments, deciding
   on a total number of bulge and flare transformations to apply,
   and filling the bstarts, bends, bpowers, and brads arrays with
   parsed or default values. nbulges is also set.

   A similar procedure is performed for the nflares, and fpowers,
   frads, and fdists arrays. */

{
  int i;
  double scratch;
  
  if (bulged->count > 0) {
    
    nbulges = intmax(5,bulged->count,
		     bulge_power->count,bulge_radius->count,
		     bulge_start->count,bulge_end->count);
    
    bstarts = (double *)(malloc(nbulges*sizeof(double)));
    bends   = (double *)(malloc(nbulges*sizeof(double)));
    bpowers = (double *)(malloc(nbulges*sizeof(double)));
    brads   = (double *)(malloc(nbulges*sizeof(double)));

    /* We now record the values for bulge starts, ends, powers
       and radii. The bulge starts and ends are normalized to lie
       in the range [0,1] with "modf". */

    for(i=0;i<nbulges;i++) {
      
      bstarts[i] = (bulge_start->count > i) ? 
	 modf(bulge_start->dval[i],&scratch):DEFAULT_BSTART;

      bends[i]   = (bulge_end->count > i) ?
	 modf(bulge_end->dval[i],&scratch):DEFAULT_BEND;

      bpowers[i]   = (bulge_power->count > i) ?
	 bulge_power->dval[i]:DEFAULT_BPOWER;

      brads[i]   = (bulge_radius->count > i) ?
	 bulge_radius->dval[i]:DEFAULT_BRAD;

    }

  }

  if (flared->count > 0) {
    
    nflares = intmax(4,flared->count,
		     flare_power->count,flare_radius->count,
		     flare_distance->count);
    
    fpowers = (double *)(malloc(nflares*sizeof(double)));
    frads   = (double *)(malloc(nflares*sizeof(double)));
    fdists  = (double *)(malloc(nflares*sizeof(double)));
 
    for(i=0;i<nflares;i++) {
      
      fdists[i] = (flare_distance->count > i) ? 
	*(flare_distance->dval):DEFAULT_FDIST;

      if (fdists[i] < 0) { fdists[i] = DEFAULT_FDIST; }
      if (fdists[i] > 1) { fdists[i] = DEFAULT_FDIST; }

      fpowers[i]   = (flare_power->count > i) ?
	*(flare_power->dval):DEFAULT_FPOWER;

      frads[i]   = (flare_radius->count > i) ?
	*(flare_radius->dval):DEFAULT_FRAD;

    }

  }

}
     

double radius_func(int comp, int vert)

/* Procedure finds the radius of the tube at the given vertex. */
/* This is complicated by the possibility of bulging and flaring */
/* in the tube, so we will have to work to find the radius. */

/* THIS FUNCTION USES THE GLOBAL "core". */

{
  double s,frac_s;
  int    i;
  double bpos,fpos,blen;
  double rad, pi = 2.0*acos(0);

  rad = radius->dval[comp % radius->count];
  
  /* We now compute the arclength position of this vertex on the tube */
  /* both in absolute terms, and as a fraction of the total arclength. */
 
  s = plc_subarc_length(core,comp,0,vert);
  frac_s = s/plc_subarc_length(core,comp,0,core->cp[comp].nv-1);
  
  /* We search the bulge intervals first. */

  for(i=0;i<nbulges;i++) {

    if (bstarts[i] < bends[i]) {  

      if (bstarts[i] <= frac_s && frac_s <= bends[i]) {

	bpos = frac_s - bstarts[i];
	bpos /= (bends[i] - bstarts[i]);

	rad += *(radius->dval)*brads[i]*  \
	  pow((sin(bpos*2.0*pi - pi/2.0) + 1)/2.0,bpowers[i]);

      }

    } else { 			/* The bulge wraps around the endpoint. */

      blen = (1.0 - bstarts[i]) + bends[i];

      if (bstarts[i] <= frac_s) {

	bpos = (frac_s - bstarts[i])/blen;
	rad += *(radius->dval)*brads[i]*  \
	  pow((sin(bpos*2.0*pi - pi/2.0) + 1)/2.0,bpowers[i]);

      } else if (frac_s <= bends[i]) {

	bpos = frac_s/blen + (1 - bstarts[i])/blen;
	rad += *(radius->dval)*brads[i]*  \
	  pow((sin(bpos*2.0*pi - pi/2.0) + 1)/2.0,bpowers[i]);

      }

    }

  }

  /* Flares come next. */
  
  for(i=0;i<nflares;i++) {

    if (core->cp[i].open) {	/* Only open components can flare. */
      
      if (frac_s < fdists[i]) {
	
	fpos = 1.0 - frac_s/fdists[i];
	rad  += *(radius->dval)*frads[i]*pow(fpos,fpowers[i]);
	
      } 
      
      if (frac_s > 1.0 - fdists[i]) {
	
	fpos = (frac_s - (1.0 - fdists[i]))/fdists[i];
	rad  += *(radius->dval)*frads[i]*pow(fpos,fpowers[i]);
	
      }
    }
  } 
  
  return rad;
}   

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  surface        *tube_surf;
  plCurve   *frameA,*frameB;

  int            infilenum,i,nerrors;
  int            *target_verts;
  plCurve   *newcore;
  plc_vector *fzbuf;

  void *argtable[] = 
    {r = arg_dbln(NULL,"cr,red","<x>",0,256,"red component of tube color"),
     g = arg_dbln(NULL,"cg,green","<x>",0,256,"green component of tube color"),
     b = arg_dbln(NULL,"cb,blue","<x>",0,256,"blue component of tube color"),
     ratio  = arg_dbl0(NULL,"stepratio","<x>","increase to add sides to tube"),
     minsides = arg_int0(NULL,"minsides","<n>","minimum number of sides for any ring of tube"),
     capped = arg_lit0("p","capped","cap the tubes with polygons"),
     closed = arg_lit0("c","closed","force tube components to close"),
     uv = arg_lit0(NULL,"uvfile","generate a u-v coordinate file for texture mapping"),
     equalized = arg_lit0("e","equalize","make core equilateral"),
     verbose = arg_lit0("v","verbose","print debugging information"),
     makerings = arg_lit0(NULL,"rings","create rings vect file (for debugging)"),
     flared  = arg_litn(NULL,"flare",0,256,"flare the ends of the tube"),
     flare_power  = arg_dbln(NULL,"fp,fpower","<x>",0,256,"controls curvature of flare"),
     flare_radius  = arg_dbln(NULL,"fr,fradius","<x>",0,256,"multiple of radius to flare to"),
     flare_distance = arg_dbln(NULL,"fd,fdist","<x>",0,256,"portion of tube to flare "
			       "(as fraction of arclength)"),
     bulged  = arg_litn(NULL,"bulge",0,256,"bulge center of tube"),
     bulge_power = arg_dbln(NULL,"bp,bpower",
			    "<x>",0,256,
			    "controls curvature of bulge"),
     bulge_radius  = arg_dbln(NULL,"br,bradius","<x>",
			      0,256,"multiple of standard radius "
			      "to bulge to"),
     bulge_start = arg_dbln(NULL,"bs,bstart","<x>",0,256,"start of bulge "
			    "(as fraction of arclength)"),
     bulge_end = arg_dbln(NULL,"be,bend","<x>",0,256,"end of bulge "
			  "(as fraction of arclength)"),
     radius = arg_dbln("r","radius","<x>",1,1000,"radius of tube component"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  fprintf(stderr,"tube (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");

  if (arg_nullcheck(argtable) != 0)
    printf("tube: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"tube");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("tube creates a Geomview OFF tube around a Geomview VECT file.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  QUIET = quiet->count > 0; /* Register if we're in batch mode */
  MAKERINGS = makerings->count > 0; /* Register if we're generating rings (for debugging) */
  
  if (minsides->count > 0) { MINSIDES = minsides->ival[0]; }
  
  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"tube: Couldn't open file %s.\n",corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }

     int plr_error_num;
     char plr_error_str[1024];
  
     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     /* We now demonstrate the octrope library's error handling protocol: */
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"tube: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);
  
     /* Preprocess the core polyline (if requested). */

     if (equalized->count > 0) {	/* Equalize sides of core polyline. */

       target_verts = (int *)(calloc(core->nc,sizeof(int)));
       for(i=0;i<core->nc;i++) target_verts[i] = core->cp[i].nv;
       newcore = plCurve_equalize_sides(core,target_verts);
       
       plc_free(core);
       core = newcore;
      
     }

     plCurve *newCore;
								       
     newCore = split_sharp_corners(core);  
     /* Preprocess by adding vertices to 
	sides incident to corners with more 
	than 5 degrees of turning angle. */
     plc_free(core);
     core = newCore;
     

     /* Force closure of components. */

     if (closed->count > 0) {	

       plc_force_closed(core);

     }

     /* Deal with bulge, flare, color arguments. */

     parse_bulge_and_flare_args();
     parse_color_args(*core);

     /* Create tube frame. */

     fzbuf = random_framezeros(core);
     plCurve_bishop_frame(core,fzbuf,&frameA,&frameB);
     free(fzbuf);

     /* Make the actual tube surface */

     if (uv->count > 0) {

       struct uvbuf uvb;

       tube_surf = make_tube(core,radius_func,
			     (ratio->count > 0 ? *(ratio->dval) : 1.0),
			     frameA,frameB,(capped->count > 0),
			     &uvb);

       FILE *uvfile;

       uvfile = fmangle(corefile->filename[infilenum],".vect",".tube.uv");
       write_uvfile(&uvb,uvfile);
       fclose(uvfile);

       if (!QUIET) {printf("tube: Wrote uv coordinate file for %d vertex, "
			   "%d face tube.\n",
			   tube_surf->verts,tube_surf->faces);}

     } else {

       tube_surf = make_tube(core,radius_func,
			     (ratio->count > 0 ? *(ratio->dval) : 1.0),
			     frameA,frameB,(capped->count > 0),
			     NULL);

     }


     if (MAKERINGS) {

       FILE *ringfile;

       ringfile = fmangle(corefile->basename[infilenum],".vect",".rings.vect");
       plc_write(ringfile,rings);
       fclose(ringfile);
       plc_free(rings);

     }

     plc_free(core);
     plc_free(frameA);
     plc_free(frameB);

     /* Post-process the tube surface. */
     /* Mangle filenames for output file. */

     char *outfile_name;

     outfile_fptr = fmangle(corefile->basename[infilenum],".vect",".tube.off");
     outfile_name = mangle(corefile->basename[infilenum],".vect",".tube.off");

     /* Now write the tube to the file. */

     write_surf_to_OFF(tube_surf,outfile_fptr);
     if (!QUIET) {printf("tube: Wrote %d vertex, %d face tube.\n",
			 tube_surf->verts,tube_surf->faces); }
     if (QUIET) {printf("%s\n",outfile_name);}

     free(outfile_name);
     fclose(outfile_fptr);
     kill_surface(tube_surf);
     free(tube_surf);

  }

  free(helpend);
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  return 0;

}




