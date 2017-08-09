/*******************************************************************

  splinevect.c : Uses splining to redistribute vertices on a plCurve
                 to obtain edges of given edgelength or according to 
                 other criteria.

	   ************************************************************/

#include<splinevect.h>
#include<argtable2.h>

#ifdef HAVE_OCTROPE

#include<octrope.h>

#endif

extern int  octrope_error_num;
extern char octrope_error_str[80];

/**********************************************************************/

struct arg_lit  *verbose;
struct arg_file *corefile;
struct arg_file *outfile;
struct arg_lit  *help;
struct arg_dbl  *edgelength;
struct arg_dbl  *resolution;
struct arg_lit  *equilateralize;
struct arg_int  *verts;
struct arg_lit  *info;
struct arg_lit  *mps;
struct arg_end  *end;
struct arg_end  *helpend;

plCurve *core;
FILE    *infile_fptr,*outfile_fptr;

/* 
 * The following procedure steps around the polygon redistributing
 * vertices as necessary to make things approach equilateral. This
 * code adapted from Eric Rawdon's TOROS polygonal runaround -- not to
 * be confused with the spherical runaround or the tangential
 * stepper. Contacts: rawdon@mathcs.duq.edu piatek@mathcs.duq.edu
 */
plCurve *plc_equilateralize( plCurve* inLink )
{
  plCurve* fixed;
  int cItr, vItr;
  double *complength = malloc(sizeof(double)*inLink->nc);
  int    *compedge = malloc(sizeof(int)*inLink->nc);
  
  fixed = plc_copy(inLink);
  plc_fix_wrap(inLink);     /* This shouldn't be needed, but... */ 

  plc_arclength(inLink,complength); /* Compute arclength for components */
  plc_edges(inLink,compedge);
  
  for( cItr=0; cItr<fixed->nc; cItr++ ) {

    double goal, used, tmpgoal, left;
    int i=0, j=1, edges;
    double* lengths;
    plc_vector* sides;
    
    edges = compedge[cItr];
    goal = complength[cItr]/(double)(edges);
        
    // we need the edge lengths and sides
    lengths = (double*)malloc(sizeof(double)*edges);
    sides = (plc_vector*)malloc(sizeof(plc_vector)*edges);
    
    for( vItr=0; vItr<edges; vItr++ ) {
      
      sides[vItr] = plc_vect_diff(inLink->cp[cItr].vt[vItr+1],
				  inLink->cp[cItr].vt[vItr]);
      lengths[vItr] = plc_M_norm(sides[vItr]);
    }
    
    left = lengths[0];
    used = 0;
    tmpgoal = goal;
    
    // commence runaround!
    while( i<edges ) {
      while( left > tmpgoal ) {
	fixed->cp[cItr].vt[j] = plc_vlincomb((tmpgoal+used)/lengths[i],sides[i],
					     1.0,inLink->cp[cItr].vt[i]);
	left -= tmpgoal;
	used += tmpgoal;
	tmpgoal = goal;
	j++;
      }
      
      tmpgoal -= left;
      i++;
      left = lengths[i];
      used = 0;
    }
    
    free(lengths);
    free(sides);
  }
  
  free(complength);
  free(compedge);

  return fixed;
}

bool resolution_ok(plCurve *inLink,int *nv) 

{
  bool ok=true;
  int i;

  for(i=0;i<inLink->nc;i++) {

    if (inLink->cp[i].nv < nv[i]) { ok = false; }
    
  }

  return ok;

}


/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,j,nerrors;
   
  void *argtable[] = 
    {
     verbose = arg_lit0(NULL,"verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     outfile = arg_filen("o","output","<file>",0,100000,"output filename(s)"),
     edgelength = arg_dbl0("e","edgelength","<x>","splines to edges of length < x"),
     resolution = arg_dbl0("r","resolution","<n>","spline to n vertices per unit length"),
     equilateralize = arg_lit0(NULL,"eq","equilateralize the splined polygon"),
     info = arg_lit0("i","info","display resolution information without splining"),
     mps = arg_lit0(NULL,"mps","minrad preserving spline (instead of cubic)"),
     verts = arg_intn("v","verts","<n>",0,100000,"number of verts (for each curve)"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};

  char plcversion[512];
  char octversion[512];
  
  plc_version(plcversion,sizeof(plcversion));

  printf("splinevect (%s, plCurve %s,",PACKAGE_STRING,plcversion);
#ifdef HAVE_OCTROPE
  octrope_version(octversion,sizeof(octversion));
  printf(" octrope version %s) \n",octversion);
#else
  printf(" without octrope)\n");
#endif
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("splinevect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"splinevect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("splinevect redistributes vertices on a VECT polyline\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (edgelength->count == 0 && resolution->count == 0 && info->count == 0 && verts->count == 0) {
  
    printf("splinevect redistributes vertices on a VECT polyline\n"
	   "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);

  }

  double res = 0.001,len = 1000;

  if (resolution->count > 0) {

    res = resolution->dval[0];

  }

  if (edgelength->count > 0) {

    len = edgelength->dval[0];

  }

  res = (res > 1.0/len) ? res : 1.0/len;

  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"splinevect: Couldn't open file %s.\n",corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }
  
     int plr_error_num;
     char plr_error_str[1024];

     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"splinevect: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     printf("splinevect: Read link from %s.\n",
	    corefile->basename[infilenum]);

     /* Now we actually do the splining... */
     
     plCurve *newCurve = NULL;
     plc_spline *spline = NULL;
     bool ok;
     double totlen,newlen,thi;
     
     double     *length = malloc(sizeof(double)*core->nc);
     int        *nv = malloc(sizeof(int)*core->nc);
     
     assert(nv != NULL);
     assert(length != NULL);
     
     totlen = plc_arclength(core,length);
     
#ifdef HAVE_OCTROPE
     thi = octrope_thickness(core,NULL,0,1.0);
     for(j=0;j<core->nc;j++) { length[j] /= thi; }
     totlen /= thi;
#else
     thi = 1.0;
#endif

     if (thi < 1e-10) { 

       fprintf(stderr,"splinevect: Thickness of curve is %g, which is too small for meaningful resolution calculation.\n",thi);
       exit(1);

     }

     /* Computes an effective resolution */

     if (verts->count > 0) {

       res = verts->ival[infilenum%verts->count]/totlen;

     }
     
     if (resolution->count != 0 || edgelength->count != 0 || verts->count != 0) {
       int totvt = 0;

       for(j=0;j<core->nc;j++) {nv[j] = ceil(res*length[j]); totvt += nv[j];}  
       
       if (verts->count > 0) { // Make things come out nicely.
	 
	 nv[j-1] += verts->ival[infilenum%verts->count] - totvt;

       }
       
       if (mps->count == 0) {
	 
	 spline = plc_convert_to_spline(core,&ok);
	 assert(ok);
	 newCurve = plc_convert_from_spline(spline,nv);
	 
       } else { 
	 
	 int failsafe;
	 plCurve *swapCurve;
	 
	 newCurve = plc_copy(core);
	 
	 for(failsafe=0;!resolution_ok(newCurve,nv);failsafe++) {
	   
	   swapCurve = plc_double_verts(newCurve);
	   
	   if (plc_num_verts(swapCurve) != 2*plc_num_verts(newCurve)) {
	     
	     printf("splinevect: mps error. Verts of doubled curve (%d) are not twice verts of old curve (%d).\n",plc_num_verts(swapCurve),plc_num_verts(newCurve));
	     exit(1);
	     
	   }
	   
	   plc_free(newCurve);
	   newCurve = plc_copy(swapCurve);
	   plc_free(swapCurve);

	   assert(failsafe < 10);
	   
	 }

       }

     } else {  /* We are getting information only. */

       newCurve = plc_copy(core);
       res = plc_num_verts(core)/totlen;

     }

     double *length_table;
     length_table = calloc(core->nc,sizeof(double));
     newlen = plc_arclength(core,length_table);
  
     printf("splinevect: Converted %g unit, %d vert plCurve ( ",totlen,plc_num_verts(core));
     for(j=0;j<core->nc;j++) { printf("%d ",core->cp[j].nv); } 
     printf(") at resolution %g",plc_num_verts(core)/totlen);

     double newthi;

#ifdef HAVE_OCTROPE
     printf(" verts/rop.\n");
     newthi = octrope_thickness(newCurve,NULL,0,1.0);
#else
     printf(" verts/len.\n");
     newthi = 1.0;
#endif

     newlen /= newthi;
     
     printf("            To %g unit, %d vert plCurve ( ",newlen,plc_num_verts(newCurve));
     for(j=0;j<newCurve->nc;j++) { printf("%d ",newCurve->cp[j].nv); } 
     printf(") at resolution %g (",plc_num_verts(newCurve)/newlen);
     for(j=0;j<newCurve->nc;j++) { printf("%g ",newCurve->cp[j].nv/(length_table[j]/newthi)); } 
     printf(")");

     free(length_table);

     char outfile_name[16000],outfile_tok[16000], *tok,*ltok;

#ifdef HAVE_OCTROPE
     printf(" verts/rop.\n");
#else
     printf(" verts/len.\n");
#endif

     /* Detect eq option and equilateralize if needed */

     if (equilateralize->count > 0) {

       double longedge,shortedge,meanedge,varedge;
       plc_edgelength_stats(newCurve,&longedge,&shortedge,&meanedge,&varedge);

       printf("max/min edgelength is %g.\n",longedge/shortedge);
       printf("equilateralizing...");

       plCurve *eqd;
       eqd = plc_equilateralize(newCurve);
       plc_edgelength_stats(eqd,&longedge,&shortedge,&meanedge,&varedge);
       printf("max/min edgelength is %g.\n",longedge/shortedge);

       if (longedge/shortedge > 1.1) { 

	 printf("equilateralization failed. quitting!.\n");
	 exit(1);

       }

       plc_free(newCurve);
       newCurve = eqd;
       
     }
     

     sprintf(outfile_tok,"%s",corefile->basename[infilenum]);

     // We need to detect the .rN token if present in filename.
     // First, we delete any trailing .vect 

     tok = strstr(outfile_tok,".vect");
     if (tok != NULL) {*tok = 0;}  // Truncate the string, killing .vect
     
     // Now search for the last . in the truncated string.

     for(ltok=NULL,tok = strstr(outfile_tok,"."); 
	 tok != NULL; 
	 ltok = tok, tok=strstr(tok+1,".")); 

     double oldres;

     if (ltok != NULL) { // We might have had only one . to start

       if (sscanf(ltok,".r%lf",&oldres) == 1) { // It was an .rN string
	 
	 *ltok = 0; // truncate it away, killing .rN

       }

     }

     // The string outfile_tok is now the truncated base of the new name

     sprintf(outfile_name,"%s.r%d.vect",outfile_tok,(int)(ceil(res)));

     // We now scrap all this if the user has overridden the automatic filename generation

     if (outfile->count > 0) {

       sprintf(outfile_name,"%s",outfile->filename[infilenum%outfile->count]);

     }

     outfile_fptr = fopen(outfile_name,"w");
       
     if (outfile_fptr == NULL) {
       
       fprintf(stderr,
	       "splinevect: Couldn't open %s for writing.\n",outfile_name);
       exit(1);
	 
     }
       
     /* Now write the component to the file. */
       
     plc_write(outfile_fptr,newCurve);
     fclose(outfile_fptr);

     printf("splinevect: Wrote link to %s.\n",
	    outfile_name);

     plc_free(newCurve);
     plc_spline_free(spline);
     free(nv);
     free(length);
     plc_free(core);    
     
  }

  if (verbose->count > 0) {

    printf("\nsplinevect processed %d input files.\n",infilenum);

  }

  return 0;

}




