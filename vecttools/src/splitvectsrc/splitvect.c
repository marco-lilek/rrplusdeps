/*******************************************************************

  splitvect.c : Divides Geomview vect files into connected components,
                or splits them in other ways.
	   
	   ************************************************************/

#include<splitvect.h>
#include<../utilib/ordie.h>
#include<../utilib/mangle.h>
#include<argtable2.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];

plCurve *plc_component(plCurve *core,int j)

/* Returns component j of curve as a new curve. */

{
  int k;
  plCurve *component;
  
  component = plc_new(1,
		      &(core->cp[j].nv),
		      &(core->cp[j].open),
		      &(core->cp[j].cc));
  
  /* We have to copy the information into component.*/
  
  for(k=0;k<component->cp[0].nv;k++) {
    
    component->cp[0].vt[k] = core->cp[j].vt[k];
    
  }
  
  for(k=0;k<component->cp[0].cc;k++) {
    
    component->cp[0].clr[k] = core->cp[j].clr[k];
    
  }

  plc_fix_wrap(component);
  
  return component;

}

plCurve *plc_split_on_long_edges(plCurve *curve,double edgelen) 

/* Divides each component of curve into several others by deleting all
   edges of length > edgelen. */

{

  bool foundlong;
  int i,j;
  int *newcmps;
  plCurve *split;

  newcmps = calloc(curve->nc,sizeof(int));

  for(i=0;i<curve->nc;i++) {

    newcmps[i]++; /* We start with one component for this one. */

    for(foundlong=false,j=0;j<curve->cp[i].nv;j++) {

      if (plc_distance(curve->cp[i].vt[j],curve->cp[i].vt[j+1]) > edgelen) {

	newcmps[i]++; /* We add another every time we delete an edge */
	foundlong = true;

      }

    }

  }

  /* We now allocate the new curve. */

  int totalcmps = 0;

  for(i=0;i<curve->nc;i++) { totalcmps += newcmps[i]; }

  int *nv, *cc;
  bool *open;

  nv = calloc(totalcmps,sizeof(int));
  cc = calloc(totalcmps,sizeof(int));
  open = calloc(totalcmps,sizeof(bool));

  int cmp=0;

  for(i=0;i<curve->nc;i++) {

    for(j=0;j<newcmps[i];j++,cmp++) {

      nv[cmp] = curve->cp[i].nv;
      open[cmp] = true;
      cc[cmp] = 0;

    }

  }

  split = plc_new(totalcmps,nv,open,cc);

  /* Now we actually fill the vertex buffers */
  
  cmp = 0; 

  for(i=0;i<totalcmps;i++) { split->cp[i].nv = 0; }

  for(i=0;i<curve->nc;i++) {

    for(j=0;j<curve->cp[i].nv;j++) {

      split->cp[cmp].vt[split->cp[cmp].nv++] = curve->cp[i].vt[j];

      if (plc_distance(curve->cp[i].vt[j],curve->cp[i].vt[j+1]) > edgelen) {

	cmp++;
	
      }

    }

  }

  plc_fix_wrap(split);
  return split;

}

int min(int a,int b) { return (a < b) ? a : b; }
int max(int a,int b) { return (a > b) ? a : b; }

int vb(plCurve *L,int cmp,int sv,int ev)

/* Number of vertices between sv and ev along curve (inclusive). This count is a little funny,
   because if ev < sv, then this refers to nothing (on an open curve) or the "back half"
   (on a closed curve). */

{
  if (sv <= ev) { 

    ev = min(L->cp[cmp].nv-1,ev);
    sv = max(0,sv);

    return ev - sv + 1;

  } else {

    if (L->cp[cmp].open) { 

      return 0;
  
    } else {

      return L->cp[cmp].nv - (sv - ev - 1);

    }

  }

}
  
plCurve *plc_slice(plCurve *L,int nfrags,int cmp[],int sv[],int ev[])

/* Creates a new plCurve with nfrags components, each containing the
   vertices from sv[i] to ev[i] (inclusive). The fragments will be
   closed if they consist of an entire component of L, which is
   denoted by setting ev > L->cp[cmp].nv-1 and sv = 0.  Constraints
   are discarded during this process in this implementation. */

{
  plCurve *newL;
  bool    *open;
  int     *cc,*nv;

  assert(nfrags > 0);
  assert(cmp != NULL && sv != NULL && ev != NULL);
  
  open = calloc_or_die(nfrags,sizeof(bool));
  cc   = calloc_or_die(nfrags,sizeof(int));
  nv   = calloc_or_die(nfrags,sizeof(int));

  int i;

  for(i=0;i<nfrags;i++) {

    if (sv[i] == 0 && ev[i] > L->cp[cmp[i]].nv-1 && !L->cp[cmp[i]].open) {

      open[i] = false;

    } else {

      open[i] = true;

    }

    if (L->cp[cmp[i]].cc == 0) { cc[i] = 0; } 
    else if (L->cp[cmp[i]].cc == 1) { cc[i] = 1; }
    else if (L->cp[cmp[i]].cc > 1) { cc[i] = vb(L,cmp[i],sv[i],ev[i]); }

    nv[i] = vb(L,cmp[i],sv[i],ev[i]);

  }

  /* Build the actual curve. */

  newL = plc_new(nfrags,nv,open,cc);

  /* Now go ahead and fill it with vertices and colors. */

  for(i=0;i<nfrags;i++) {

    int vt,run;
    run = vb(L,cmp[i],sv[i],ev[i]);

    for(vt=sv[i];vt<sv[i]+run;vt++) {

      newL->cp[i].vt[vt-sv[i]] = L->cp[cmp[i]].vt[(vt % L->cp[cmp[i]].nv)];
      
    }

    if (newL->cp[i].cc == 1) { 

      newL->cp[i].clr[0] = L->cp[cmp[i]].clr[0];

    } else if (newL->cp[i].cc > 1) {

      for(vt=sv[i];vt<sv[i]+run;vt++) {

	newL->cp[i].clr[vt-sv[i]] = L->cp[cmp[i]].clr[(vt % L->cp[cmp[i]].nv)];
      
      }

    }

  }

  return newL;
  
}

/**********************************************************************/

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_dbl  *splitlong;
struct arg_end *end;
struct arg_end *helpend;
struct arg_str  *vertrange;

plCurve *core;
FILE         *infile_fptr,*outfile_fptr;

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,j,nerrors;
  plCurve       *component;
  char           outfile_name[1000],outfile_tail[1000];
 
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     splitlong = arg_dbl0(NULL,"split-on-edgelength","<x>","splits at edges of length > x"),
     vertrange = arg_strn(NULL,"verts","c-d,n",0,100000,"extract verts c-d (inclusive) of component n. c,d can be 'first' or 'last'"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("splitvect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"splitvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("splitvect divides a VECT file into its component polylines\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  /* For slicing purposes, we now parse the vertex range arguments */

  int *cmp,*sv,*ev;

  cmp = calloc_or_die(vertrange->count,sizeof(int));
  sv  = calloc_or_die(vertrange->count,sizeof(int));
  ev  = calloc_or_die(vertrange->count,sizeof(int));

 
  
  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"splitvect: Couldn't open file %s.\n",corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }
  
     int plr_error_num;
     char plr_error_str[1024];

     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"splitvect: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     /* Now we decide what to do, based on the arguments. */

     if (vertrange->count > 0) {
       
       int i;
       
       for (i=0;i<vertrange->count;i++) { /* For all ranges to constrain. */
	 
	 /* Start by parsing the range. */
	 
	 if (sscanf(vertrange->sval[i],"%d-%d,%d",
		    &(sv[i]),&(ev[i]),&(cmp[i])) == 3) {
	   
	   assert(cmp[i] < core->nc);
	   
	   /* Do nothing, we're good. */
	   
	 } else if (sscanf(vertrange->sval[i],"first-%d,%d",
			   &(ev[i]),&(cmp[i])) == 2) {
	   
	   assert(cmp[i] < core->nc);
	   sv[i] = 0;
	   
	 } else if (sscanf(vertrange->sval[i],"%d-last,%d",
			   &(sv[i]),&(cmp[i])) == 2) {
	   
	   assert(cmp[i] < core->nc);
	   ev[i] = core->cp[cmp[i]].nv-1;
	   
	 } else if (sscanf(vertrange->sval[i],"first-first,%d",
			   &(cmp[i])) == 1) {
	   
	   assert(cmp[i] < core->nc);
	   sv[i] = ev[i] = 0;
	   
	 } else if (sscanf(vertrange->sval[i],"last-last,%d",
			   &(cmp[i])) == 1) {
	   
	   assert(cmp[i] < core->nc);
	   sv[i] = ev[i] = core->cp[cmp[i]].nv-1;
	   
	 } else if (sscanf(vertrange->sval[i],"first-last,%d",
			   &(cmp[i])) == 1) {
	   
	   assert(cmp[i] < core->nc);
	   sv[i] = 0; ev[i] = core->cp[cmp[i]].nv-1;
	   
	 } else {
	   
	   printf("constrainvect: Couldn't parse %s.\n",
		  vertrange->sval[i]);
	   
	   exit(1);
	   
	 }
	 
	 if (verbose->count > 0) {
	   
	   printf("constrainvect: Parsed vertrange %s to mean %d verts from (%d,%d).\n",
		  vertrange->sval[i],vb(core,cmp[i],sv[i],ev[i]),sv[i],ev[i]);
	   
	 }
	 
       }
       
       FILE *outfile;
       char *outname;

       outfile = fmangle(corefile->basename[infilenum],"vect","range.vect");
       outname = mangle(corefile->basename[infilenum],"vect","range.vect");

       plCurve *slice;

       slice = plc_slice(core,vertrange->count,cmp,sv,ev);

       plc_write(outfile,slice);
       fclose(outfile);
       plc_free(slice);

       printf("splitvect: Wrote output file to %s.\n",outname);
       free(outname);

     } else {

       /* Now we divide the link into components. */
       
       plCurve *split;
       
       if (splitlong->count > 0) {
       
	 printf("splitvect: Now deleting edges longer than %g.\n",splitlong->dval[0]);
	 split = plc_split_on_long_edges(core,splitlong->dval[0]);
	 plc_free(core);
	 core = split;
	 
       } 
       
       for(j=0;j<core->nc;j++) {
	 
	 component = plc_component(core,j);
	 
       /* Now we write the component to a file. */
	 
	 if (strlen(corefile->basename[infilenum]) > sizeof(outfile_name)-20) {
	   
	   fprintf(stderr,"splitvect: Ridiculously long input filename "
		   "can't be parsed.\n");
	   exit(1);
       
	 }
	 
	 sprintf(outfile_name,"%s",corefile->basename[infilenum]);
	 
	 if (strstr(outfile_name,".vect") != NULL) {
	   
	   sprintf(strstr(outfile_name,".vect"),".%d.vect",j+1);
	   
	 } else {
	   
	   sprintf(outfile_tail,".%d.vect",j+1);	 
	   strcat(outfile_name,outfile_tail);
	   
	 }
	 
	 outfile_fptr = fopen(outfile_name,"w");
	 
	 if (outfile_fptr == NULL) {
	   
	   fprintf(stderr,
		   "splitvect: Couldn't open %s for writing.\n",outfile_name);
	   exit(1);
	   
	 }
	 
	 /* Now write the component to the file. */
	 
	 plc_write(outfile_fptr,component);
	 fclose(outfile_fptr);
	 plc_free(component);
	 
       }
       
       plc_free(core);
       
       printf("splitvect: Produced %d component files from %s.\n",
	      j,corefile->basename[infilenum]);

     }
     
  }

  if (verbose->count > 0) {

    printf("\nsplitvect processed %d input files.\n",infilenum);

  }

  return 0;

}




