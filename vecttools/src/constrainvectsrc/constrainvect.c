/*******************************************************************

 constrainvect.c : Moves or rotates Geomview VECT files in space.
	   
	   ************************************************************/

#include<constrainvect.h>
#include<argtable2.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];

/**********************************************************************/

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_lit  *endpoints;
struct arg_lit  *help;
struct arg_str  *planes;
struct arg_lit  *fix;
struct arg_lit  *show;
struct arg_lit  *force;
struct arg_lit  *geomview;
struct arg_lit  *clear;
struct arg_str  *vertrange;
struct arg_end *end;
struct arg_end *helpend;

plCurve *core;
FILE         *infile_fptr,*outfile_fptr;

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,nerrors;
  int            i;
 
  void *argtable[] = 
    {verbose = arg_lit0(NULL,"verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     vertrange = arg_strn("v","verts","c-d,n",0,100000,"constrain verts c-d (inclusive) of component n. c,d can be 'first' or 'last'"),
     fix = arg_lit0("f","fix","fix at current location (only applies to first vert in range)"),
     endpoints = arg_lit0(NULL,"endpoints","apply constraint to (all) endpoints"),
     planes = arg_strn("p","plane","x,y,z",0,100000,"fix to plane with normal (x,y,z) containing first vert in range"),
     show = arg_lit0("s","show","show constraints"),
     force = arg_lit0(NULL,"force","force curve to obey constraints"),
     geomview = arg_lit0("g","geomview","draw constraints in geomview"),
     clear = arg_lit0("c","clear","remove all constraints"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("constrainvect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"constrainvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("constrainvect adds constraints to a VECT curve\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }
  
  if (fix->count > 0 && planes->count > 0) {

    printf("constrainvect: Can't give both fixed and plane constraints\n"
	   "               in the same run. You must run several times in\n"
	   "               order to mix constraints (it is ok to run on\n"
	   "               the same file multiple times).\n\n");
    
    exit(0);

  }

  /* Now we parse the plane constraints, if they exist. */

  plc_vector *normals;
  normals = calloc(planes->count,sizeof(plc_vector));

  if (planes->count > 0) {

    assert(normals!=NULL); 

    for(i=0;i<planes->count;i++) {

      if (sscanf(planes->sval[0],"%lf,%lf,%lf",
		 &(normals[i].c[0]),&(normals[i].c[1]),&(normals[i].c[2])) != 3) {
	
	printf("constrainvect: Couldn't parse %s.\n",
	       planes->sval[0]);

	exit(1);
	
      }

      bool ok;
      
      normals[i] = plc_normalize_vect(normals[i],&ok);

      if (!ok) { 

	printf("constrainvect: Normal to plane must be nonzero.\n");
	exit(1);

      }
 
    }

  }
   
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

     if (show->count > 0) { 

       plc_constraint *this;

       for(i=0,this=core->cst;this!=NULL;this=this->next,i++); 
       /* Count constraints */

       printf("constrainvect: File %s has %d constraints.\n",
	      corefile->basename[infilenum],i);

       for(this=core->cst;this!=NULL;this=this->next) {

	 printf("  Component: %2d Verts: %3d-%3d Type: ",
		this->cmp,this->vert,this->vert+this->num_verts-1);
	 
	 if (this->kind == unconstrained) {
	   
	   printf("UNCONSTRAINED\n");

	 } else if (this->kind == fixed) {

	   printf("FIXED\n");

	 } else if (this->kind == line) {

	   printf("LINE Dir: (%4g,%4g,%4g) Point: (%4g,%4g,%4g)\n",
		  plc_M_clist( this->vect[0] ), plc_M_clist( this->vect[1] ));

	 } else if (this->kind == plane) {

	   printf("PLANE Normal: (%4g,%4g,%4g) Dist: %4g\n",
		  plc_M_clist( this->vect[0] ), this->vect[1].c[0]);

	 }

       }

     }

     if (geomview->count > 0) { 

       plc_constraint *this;
       plc_constraint points[1024],planes[1024],lines[1024];
       int            npoints=0,nplanes=0,nlines=0;
       
       for(i=0,this=core->cst;this!=NULL;this=this->next,i++); 
       /* Count constraints */

       printf("constrainvect: File %s has %d constraints.\n",
	      corefile->basename[infilenum],i);

       for(this=core->cst;this!=NULL;this=this->next) {

	 if (this->kind == unconstrained) {
	   
	   // Do nothing.

	 } else if (this->kind == fixed) {

	   points[npoints++] = *this;

	 } else if (this->kind == line) {

	   lines[nlines++] = *this;

	 } else if (this->kind == plane) {

	   planes[nplanes] = *this;
	   /* Now we modify this to put the plane's vect array into normal, point form. */
	   bool ok;

	   planes[nplanes].vect[0] = plc_normalize_vect(this->vect[0],&ok); assert(ok);
	   planes[nplanes].vect[1] = plc_scale_vect(this->vect[1].c[0],planes[nplanes].vect[0]);

	   nplanes++;
	   
	 }

       }

       /* We have now sorted the constraints. We start a file to display them. */
       
       char outfile_name[1024],outfile_tail[1024];

       if (strlen(corefile->basename[infilenum]) > sizeof(outfile_name)-20) {
       
	 fprintf(stderr,"constrainvect: Ridiculously long input filename "
		 "can't be parsed.\n");
	 exit(1);
       
       }
     
       sprintf(outfile_name,"%s",corefile->basename[infilenum]);
     
       if (strstr(outfile_name,".vect") != NULL) {
       
	 sprintf(strstr(outfile_name,".vect"),".list");
	 
       } else {
	 
	 sprintf(outfile_tail,".list");	 
	 strcat(outfile_name,outfile_tail);
	 
       }

       outfile_fptr = fopen(outfile_name,"w");
     
       if (outfile_fptr == NULL) {
	 
	 fprintf(stderr,"constrainvect: Couldn't open file %s.\n",
		 corefile->filename[infilenum]);
	 continue;  /* Try the next file */
	 
       }
       
       /* We now display things in order. For the points, this is easy.*/

       fprintf(outfile_fptr,"LIST\n\n");

       int i;

       for(i=0;i<npoints;i++) {

	 fprintf(outfile_fptr,"{ SPHERE 0.1 %g %g %g }\n",plc_M_clist( points[i].vect[0] ));

       }

       printf("constrainvect: Displayed %d fixed point constraints as spheres.\n",npoints);

       plCurve *planeSkel;
       int nv = 1,cc = 0;
       bool open = true;

       planeSkel = plc_new(1,&nv,&open,&cc);

       for(i=0;i<nplanes;i++) {

	 /* For each plane, try to find intersections with other planes. */

	 plc_vector intersections[1024];
	 bool ok;
	 int nintersections=0;

	 int j,k;

	 for(j=0;j<nplanes;j++) {

	   for(k=0;k<nplanes;k++) {

	     if ((i != j) && (j != k) && (i != k)) {

	       intersections[nintersections] 
		 = plc_3plane_intersection( planes[i].vect[0], planes[i].vect[1],
					    planes[j].vect[0], planes[j].vect[1],
					    planes[k].vect[0], planes[k].vect[1], &ok);

	       if (ok) { nintersections++; }

	     }

	   }

	 }

	 /* Now add them to our skeleton plCurve. */

	 bool open=false;
	 int cc=0;
	 
	 plc_add_component(planeSkel,1,nintersections,open,cc,intersections,NULL);

       }

       plc_drop_component(planeSkel,0); /* Drop the first, placeholder, component. */

       if (planeSkel->nc > 0) {

	 printf("constrainvect: Displayed %d plane constraints by intersection points.\n",planeSkel->nc);

	 fprintf(outfile_fptr,"\n{\n");
	 plc_write(outfile_fptr,planeSkel);
	 fprintf(outfile_fptr,"\n}\n");
	 plc_free(planeSkel);

       }

       /* Now we ought to do something with lines, so we do. */

       for(i=0;i<nlines;i++) {

	 printf("Warning: Did not display line constraint.\n");

       }

       fclose(outfile_fptr);

       printf("constrainvect: Geomview picture of constraints written to %s.\n",outfile_name);
       
       
     }



     /* We have showed constraints. If we need another action, do it. */

     if (clear->count > 0) {

       plc_remove_all_constraints(core);

       outfile_fptr = fopen(corefile->filename[infilenum],"w");
       
       if (outfile_fptr == NULL) {
    
	 fprintf(stderr,"constrainvect: Couldn't open file %s.\n",
		 corefile->filename[infilenum]);
	 continue;  /* Try the next file */
	 
       }
       

       plc_write(outfile_fptr,core);

       printf("constrainvect: CLEARED all constraints from %s and saved it.\n",
	      corefile->basename[infilenum]);

       continue; /* On to the next file */

     }

     for (i=0;i<vertrange->count;i++) { /* For all ranges to constrain. */

       /* Start by parsing the range. */

       int cmp;
       int start,end,num;

       if (sscanf(vertrange->sval[i],"%d-%d,%d",
		 &(start),&(end),&(cmp)) == 3) {
	
	 /* Do nothing, we're good. */

       } else if (sscanf(vertrange->sval[i],"first-%d,%d",
		  &(end),&(cmp)) == 2) {
	
	 start = 0;
	
       } else if (sscanf(vertrange->sval[i],"%d-last,%d",
		  &(start),&(cmp)) == 2) {
	
	 end = core->cp[cmp].nv-1;
	
       } else if (sscanf(vertrange->sval[i],"first-first,%d",
		  &(cmp)) == 1) {
	
	 start = end = 0;
	
       } else if (sscanf(vertrange->sval[i],"last-last,%d",
		 &(cmp)) == 1) {
	
	 start = end = core->cp[cmp].nv-1;
	
       } else if (sscanf(vertrange->sval[i],"first-last,%d",
		  &(cmp)) == 1) {
	
	 start = 0; end = core->cp[cmp].nv-1;
	
       } else {

	 printf("constrainvect: Couldn't parse %s.\n",
		vertrange->sval[i]);
	 
	 exit(1);

       }

       num = end - start + 1;

       if (num <= 0) {

	 fprintf(stderr,"constrainvect: Vertex range %s must go in order.\n",
		 vertrange->sval[i]);

       }

       if (verbose->count > 0) {

	 printf("constrainvect: Parsed vertrange %s to mean %d verts from (%d,%d).\n",
		vertrange->sval[i],num,start,cmp);

       }

       /* We have read the vertex range. Now apply it.*/ 

       if (fix->count > 0) {

	 plc_set_fixed(core,cmp,start,core->cp[cmp].vt[start]);

	 if (verbose->count > 0) {

	   printf("constrainvect: Set this range to FIXED.\n");

	 }

       } else if (planes->count > 0) {

	 plc_constrain_to_plane(core,cmp,start,num,
				normals[i % planes->count],
				plc_dot_prod(normals[i % planes->count],
					     core->cp[cmp].vt[start]));

	 if (verbose->count > 0) {

	   printf("constrainvect: Set this range to PLANE constraint %d.\n",
		  i%planes->count);

	 }
	 
       }

     }

     if (endpoints->count > 0) {

       if (fix->count > 0) { 

	 int cp;

	 for(cp=0;cp<core->nc;cp++) {

	   plc_set_fixed(core,cp,0,core->cp[cp].vt[0]);
	   plc_set_fixed(core,cp,core->cp[cp].nv-1,core->cp[cp].vt[core->cp[cp].nv-1]);

	 }

	 vertrange->count = core->nc*2;

       } else {

	 printf("--endpoints only works with --fix constraints\n");
	 exit(1);

       }

     }

     /* Now force the curve to obey the constraints, if needed */

     if (force->count > 0) {

       plc_fix_cst(core);
       printf("constrainvect: Forced curve to obey all constraints.\n");

     }

     /* We have applied all constraints. Save the file. */

     printf("constrainvect: %d constraints applied to %s.\n",
	    vertrange->count,corefile->basename[infilenum]);
     
     outfile_fptr = fopen(corefile->filename[infilenum],"w");
     
     if (outfile_fptr == NULL) {
       
       fprintf(stderr,"constrainvect: Couldn't open file %s.\n",
	       corefile->filename[infilenum]);
       continue;  /* Try the next file */
       
     }
       
     plc_write(outfile_fptr,core);

     printf("constrainvect: Saved file to %s.\n",
	    corefile->basename[infilenum]);
     

     plc_free(core);
     
     if (normals != NULL) { free(normals); }

  }

  return 0;

}




