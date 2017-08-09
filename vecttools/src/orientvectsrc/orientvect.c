/*******************************************************************

orientvect.c : Rotates plCurve to new orientation (either with respect to original axes or the inertial axes).
	   
	   ********************************************************/

#include<orientvect.h>
#include<argtable2.h>
#include<../utilib/mangle.h>
#include<../utilib/ordie.h>

extern int  octrope_error_num;
extern char octrope_error_str[80];

/******************************************************************/

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_int  *arg_inertial_axis;
struct arg_lit  *to_straights;
struct arg_str  *arg_matrix;
struct arg_dbl  *arg_spin;
struct arg_file *outfile;
struct arg_end *end;
struct arg_end *helpend;

plCurve  *core;
FILE     *infile_fptr,*outfile_fptr;
int       axis = 1;

/******************************************************************/

void inertial_axes_and_com(plCurve *L,plc_vector *center_of_mass,plc_vector *inertial_axes)

{
  int cp;
  int vt;
  int nv;
  int i;
  int j;

  double *(cloud[3]);

  /* Construct point cloud */
 
  nv = plc_num_verts(L);
  cloud[0] = calloc_or_die(nv,sizeof(double));
  cloud[1] = calloc_or_die(nv,sizeof(double));
  cloud[2] = calloc_or_die(nv,sizeof(double));

  for(cp=0,i=0;cp<L->nc;cp++) {
    for(vt=0;vt<L->cp[cp].nv;vt++,i++) {
      for(j=0;j<3;j++) {
	(cloud[j])[i] = L->cp[cp].vt[vt].c[j];
      }
    }
  }

  /* Now construct center of mass */
  
  (*center_of_mass) = plc_build_vect(0,0,0);
  for(i=0;i<nv;i++) {
    for(j=0;j<3;j++) {
      center_of_mass->c[j] += (cloud[j])[i];
    }
  }
  *center_of_mass = plc_scale_vect(1.0/nv,*center_of_mass);

  /* Now subtract that from the point cloud */

  for(i=0;i<nv;i++) {
    for(j=0;j<3;j++) {
      (cloud[j])[i] -= center_of_mass->c[j];
    }
  }
    
  /* Now construct the Gram matrix */

  double gram[3][3];
  int k;

  for(i=0;i<3;i++) {
    for(j=0;j<=i;j++) {
      gram[i][j] = 0;
      for(k=0;k<nv;k++) {
	gram[i][j] += (cloud[i])[k] * (cloud[j])[k];
      }
    }
  }

  for(i=0;i<3;i++) {
    for(j=i+1;j<3;j++) {
      gram[i][j] = gram[j][i];
    }
  }

  /* We now need to find the eigenvalues of the Gram matrix */

  double data[9];

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      data[3*i+j] = gram[i][j];
    }
  }
     
  gsl_matrix_view m 
    = gsl_matrix_view_array (data, 3, 3);
  
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_matrix *evec = gsl_matrix_alloc (3, 3);
     
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (3);
       
  gsl_eigen_symmv (&m.matrix, eval, evec, w);     
  gsl_eigen_symmv_free (w);    
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  /* We now read off the eigenvectors */

  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      inertial_axes[j].c[i] = gsl_matrix_get(evec,i,j);
    }
  }

  /* Now clean up after ourselves */ 

  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  free(cloud[0]);
  free(cloud[1]);
  free(cloud[2]);

}

void plc_orient_to_inertial_axes(plCurve *L,int axis)
{
  plc_vector inertial_axes[3];
  plc_vector center_of_mass;
  int axisorder[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
  int j;

  inertial_axes_and_com(L,&center_of_mass,inertial_axes);

  printf("Found center of mass and inertial axes:\n\n");

  printf("Center of mass: (%g,%g,%g)\n",plc_M_clist(center_of_mass));
  printf("Axis 1:      (%10g,\t%10g,\t%10g)\n",plc_M_clist(inertial_axes[0]));
  printf("Axis 2:      (%10g,\t%10g,\t%10g)\n",plc_M_clist(inertial_axes[1]));
  printf("Axis 3:      (%10g,\t%10g,\t%10g)\n",plc_M_clist(inertial_axes[2]));

  printf("Now orienting curve to inertial axes...");

  /* First, we translate the curve so the center of mass is at origin. */

  int cp, vt;

  for(cp=0;cp<L->nc;cp++) {
    for(vt=0;vt<L->cp[cp].nv;vt++) {
      plc_M_sub_vect(L->cp[cp].vt[vt],center_of_mass);
    }
  }
      
  /* Now we rotate to the new axes */

  plc_vector newvt;

  /* First, we force the frame to be right-handed */

  if (plc_dot_prod(inertial_axes[0],plc_cross_prod(inertial_axes[1],inertial_axes[2])) < 0) {

    inertial_axes[axisorder[axis][2]] = plc_scale_vect(-1,inertial_axes[axisorder[axis][2]]);

  }

  for(cp=0;cp<L->nc;cp++) {
    for(vt=0;vt<L->cp[cp].nv;vt++) {
      for(j=0;j<3;j++) {
	newvt.c[j] = plc_dot_prod(inertial_axes[axisorder[axis][j]],L->cp[cp].vt[vt]);
      }
      L->cp[cp].vt[vt] = newvt;
    }
  }

  /* Now we translate the curve back to the original center of mass. */

  for(cp=0;cp<L->nc;cp++) {
    for(vt=0;vt<L->cp[cp].nv;vt++) {
      plc_M_add_vect(L->cp[cp].vt[vt],center_of_mass);
    }
  }

  plc_fix_wrap(L);
 
  printf("done\n");
}

void plc_apply_euler_angles(plCurve *L, double spin[3])

// Applies transformation by Euler angles (in order around X, Y, Z axes) to plCurve.

{
  plc_vector workingvect,newvt;
  int cp,vt;

  printf("Now applying Euler angles %g:%g:%g ...",spin[0],spin[1],spin[2]);

  for(cp=0;cp<L->nc;cp++) {
    for(vt=0;vt<L->cp[cp].nv;vt++) {

      newvt = L->cp[cp].vt[vt];
      
      workingvect.c[1] = cos(spin[0]) * newvt.c[1] + sin(spin[0]) * newvt.c[2];
      workingvect.c[2] = -sin(spin[0]) * newvt.c[1] + cos(spin[0]) * newvt.c[2];
      workingvect.c[0] = newvt.c[0];
      newvt = workingvect;
      
      workingvect.c[2] = cos(spin[1]) * newvt.c[2] + sin(spin[1]) * newvt.c[0];
      workingvect.c[0] = -sin(spin[1]) * newvt.c[2] + cos(spin[1]) * newvt.c[0];
      workingvect.c[1] = newvt.c[1];
      newvt = workingvect;

      workingvect.c[0] = cos(spin[2]) * newvt.c[0] + sin(spin[2]) * newvt.c[1];
      workingvect.c[1] = -sin(spin[2]) * newvt.c[0] + cos(spin[2]) * newvt.c[1];
      workingvect.c[2] = newvt.c[2];
      newvt = workingvect;
      
      L->cp[cp].vt[vt] = newvt;
    }
  }
  
  plc_fix_wrap(L);

  printf("done\n");
}

void plc_apply_matrix(plCurve *L,double M[3][3])

{
  plc_vector workingvect,newvt;
  int cp,vt;
  int i,j;

  printf("Now applying matrix \n"
	 "\t %3.4g \t %3.4g \t %3.4g \n"
	 "\t %3.4g \t %3.4g \t %3.4g \n"
	 "\t %3.4g \t %3.4g \t %3.4g ... ",
	 (M[0][0]), (M[0][1]), (M[0][2]), 
	 (M[1][0]), (M[1][1]), (M[1][2]), 
	 (M[2][0]), (M[2][1]), (M[2][2])
	 );

  for(cp=0;cp<L->nc;cp++) {
    for(vt=0;vt<L->cp[cp].nv;vt++) {

      newvt = L->cp[cp].vt[vt];
      workingvect.c[0] = workingvect.c[1] = workingvect.c[2] = 0;
 
     for(i=0;i<3;i++) {
	for(j=0;j<3;j++) {
	  workingvect.c[i] += newvt.c[j]*M[i][j];
	}
     }
      
     L->cp[cp].vt[vt] = workingvect;
    }
  }
  
  plc_fix_wrap(L);

  printf("done\n");
}
  


/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,nerrors;
  double         spin[3] = {0,0,0};
  double         M[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
   
  void *argtable[] = 
    {verbose = arg_lit0(NULL,"verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     arg_inertial_axis = arg_int0("a","axis","<n>","axis of inertia (1, 2, or 3)"),
     arg_spin = arg_dbln("s","spin","<rad>",0,3,"angle of spin around axis n (may be repeated)"),
     arg_matrix = arg_str0("m","matrix","{{a11,a12,a13},{....},{a31,a32,a33}}","orientation matrix"),
     outfile = arg_filen("o","outfile","<file>",0,100000,"output files"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("orientvect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"orientvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("orientvect orients a plCurve to its axes of inertia\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (arg_inertial_axis->count > 0) { axis = arg_inertial_axis->ival[0]; }
  if (arg_spin->count > 0) { int i; for(i=0;i<arg_spin->count;i++) {spin[i%3] = arg_spin->dval[i];} }
  if (arg_matrix->count > 0) {
    if (sscanf(arg_matrix->sval[0],"{{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}}",
	       &(M[0][0]), &(M[0][1]), &(M[0][2]), 
	       &(M[1][0]), &(M[1][1]), &(M[1][2]), 
	       &(M[2][0]), &(M[2][1]), &(M[2][2])) != 9) {
      printf("orientvect: Couldn't parse '--matrix' argument of %s\n",arg_matrix->sval[0]);
      exit(1);
    }
  }
  
  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"orientvect: Couldn't open file %s.\n",
	       corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }
  
     int plr_error_num;
     char plr_error_str[1024];

     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"orientvect: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     /* Now we work on the link. */

     if (arg_inertial_axis->count > 0) {
       plc_orient_to_inertial_axes(core,axis-1);
     }

     if (arg_spin->count > 0) {
       plc_apply_euler_angles(core,spin);
     }

     if (arg_matrix->count > 0) {
       plc_apply_matrix(core,M);
     }
       
     /* We have applied all transformations. Save the file. */

     printf("orientvect: oriented %s to axis %d\n",
	    corefile->basename[infilenum],axis);
     
     outfile_fptr = NULL;

     if (outfile->count > infilenum) {

       outfile_fptr = fopen(outfile->basename[infilenum],"w");

       if (outfile_fptr != NULL) {

	 printf("orientvect: Saved file to %s.\n",outfile->basename[infilenum]);

       }
       
     } 

     if (outfile_fptr == NULL) {

       char *newname;
       
       newname = mangle(corefile->filename[infilenum],
			".vect",".oriented.vect");
       
       outfile_fptr = fmangle(corefile->filename[infilenum],
			      ".vect",".oriented.vect");

       printf("orientvect: Saved file to %s.\n",newname);

       free(newname);

     }
       
     plc_write(outfile_fptr,core);    
     plc_free(core);
     
  }

  return 0;

}




