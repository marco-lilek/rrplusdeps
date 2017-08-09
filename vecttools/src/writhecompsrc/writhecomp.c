/**************************************************************************/
/*		      							  */
/* writhecomp :This program is designed to compute the writhing number    */
/*	       of polygonal curves specified as Geomview polyline objects */
/*	       in the .VECT format using arbitrary precision arithmetic.  */
/*                                                                        */
/*             The essential problem of writhecomp is that the exact      */
/*             Banchoff algorithm for the writhing number of a polygonal  */
/*             is not particularly well-behaved numerically.              */
/*                                                                        */
/*             Note that other algorithms may indeed perform better.      */
/*             For example, the linking number based algorithm of Klein   */
/*             and Langowski seems like it presents different numerical   */
/*             challenges.                                                */
/*                                                                        */
/*             Thus, to combat the problems of bad trig functions and     */
/*             roundoff error, the program is capable of computing writhe */
/*             in three arithmetic models:                                */
/*                                                                        */
/*                   C double mode  (when speed is most important)        */
/*                   MAPM mode (arbitrary precision floating point)       */
/*                   PARI mode (fast arbitrary precision/hard to support) */
/*                                                                        */
/*             One can be confident in the results of the computation     */
/*             if they are not affected by a change in precision.         */
/*             One can also run the test suite (verifywc) to make sure    */
/*             of the accuracy of one's computations.                     */
/*                                                                        */
/*             The core of MAPM mode is distributed with writhecomp, PARI */
/*             mode will require one to compile and install GP-PARI, which*/
/*             is often in and of itself an adventure, depending on your  */
/*             architecture.                                              */
/*                                                                        */
/*             We also require:                                           */
/*                                                                        */
/*                   argtable                                             */
/*                                                                        */
/*             and we include the cubature code of Steven Johnson.        */
/*                                                                        */
/*           As of 11/2009, writhecomp has been ported to the plCurve API.*/
/*                                                                        */
/*                                                                        */
/**************************************************************************/

#include"writhecomp.h"

double fast_writhe(plCurve *L); /* An external function in cdcuhre_writhe.c */
double fast_spline_writhe(plCurve *L);

int TIMING = {FALSE};   /* Just estimate the time for a computation? no. */
int SAFE = {FALSE};	        /* Safe mode self-checking?  Default- no.*/ 
int VERBOSE = {FALSE};	/* Report everything to the user? Default- no. */
int DEC_PRECISION = {16};   /* Digits of (decimal) precision to use. */

/*************************************************************/

double estimated_time(int verts)

    /* Procedure returns the estimated time (in seconds) that it
       will take to compute the writhe of a polyline with <verts> vertices. */
     
{
  time_t t_start;
  time_t t_end;
  
  int    n,N_ITERATIONS = 1000;
  int    i;
  
  NUMVECTOR   rand_vector[4];
  NUMPTR      tmp;
  
  /* First, set up some fake vectors */
  
  rand_vector[0][0] = (*double_to_np)(0);
  rand_vector[0][1] = (*double_to_np)(0);
  rand_vector[0][2] = (*double_to_np)(0);
  
  rand_vector[1][0] = (*double_to_np)(1);
  rand_vector[1][1] = (*double_to_np)(0);
  rand_vector[1][2] = (*double_to_np)(0);
  
  rand_vector[2][0] = (*double_to_np)(0);
  rand_vector[2][1] = (*double_to_np)(1);
  rand_vector[2][2] = (*double_to_np)(0);
  
  rand_vector[3][0] = (*double_to_np)(0);
  rand_vector[3][1] = (*double_to_np)(0);
  rand_vector[3][2] = (*double_to_np)(1);
  
  /* Now compute some dihedral sums */
  
  time(&t_start);
  
  for(n=0;n<N_ITERATIONS;n++) {
    
    tmp = (*np_dihedral_sum)(rand_vector[0],rand_vector[1],
			     rand_vector[2],rand_vector[3]);
    free(tmp);
    
  }
  
  time(&t_end);
  
  /* Free the fake vectors (even small memory leaks are bad!). */
    
  for(i=0;i<4;i++) {
    
    kill_nv(rand_vector[i]);
    
  }
  
  /* Now extrapolate to guess the total time for the computation. */
  
  return (double)( ((int)difftime(t_end,t_start))*(verts*verts - verts)/(N_ITERATIONS) );
  
}

struct arg_lit  *arg_verbose;
struct arg_file *arg_infile;
struct arg_lit  *arg_safe;
struct arg_lit  *arg_fast;
struct arg_lit  *arg_fast_spline;
struct arg_str  *arg_model;
struct arg_int  *arg_precision;
struct arg_dbl  *arg_resolution;
struct arg_lit  *arg_timing;
struct arg_lit  *help;
struct arg_end  *end;
struct arg_end  *helpend;

int main(int argc, char *argv[])
{
  int nerrors;
  enum MODEL {cdouble,mapm} arithmodel = {cdouble};
  plCurve *writhe_pline;
  FILE *infile;
  
  double pwrithe = 0, lpwrithe = 0, hpwrithe = 0;

  void *argtable[] = {

    arg_infile   = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
    arg_verbose  = arg_lit0(NULL,"verbose","print debugging information"),
    arg_model    = arg_str0("m","model","<double or mapm>","arithmetic model for calculation"),
    arg_fast     = arg_lit0("f","fast","compute using fast subdivision model"),
    arg_fast_spline = arg_lit0(NULL,"fastspline","compute using (splined) fast subdivision model"),
    arg_resolution = arg_dbl0("r","resolution","<r>","perform computation at resolution r verts/unit length"),
    arg_safe     = arg_lit0(NULL,"safe","safe mode calculates in several different precisions"),
    arg_precision = arg_int0("p","precision","<digits>","(decimal) digits of precision"),
    arg_timing   = arg_lit0(NULL,"timing","estimate time for calculation (but do not compute)"),
    help         = arg_lit0("?","help","display help message"),
    end          = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
 
    /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("writhe: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"writhe");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("writhe computes the writhing number of a VECT polyline\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  /* We've parsed the arguments. Now slot them into the proper globals. */

  VERBOSE = (arg_verbose->count > 0);
  
  if (arg_model->count > 0) {

    if (!strcmp(arg_model->sval[0],"double")) { arithmodel = cdouble; }
    else if (!strcmp(arg_model->sval[0],"mapm")) { arithmodel = mapm; }

  }

  SAFE = (arg_safe->count > 0); 

  if (SAFE) { arithmodel = mapm; }

  if (arg_precision->count > 0) {

    DEC_PRECISION = arg_precision->ival[0];

  } else {

    DEC_PRECISION = 16;

  }

  TIMING = (arg_timing->count > 0); 
    
  // We have processed the command-line arguments and must act on them.
  // Our first task is to see whether they make sense.
  
  int infilenum;

  for(infilenum=0;infilenum < arg_infile->count;infilenum++) {

    infile = fopen(arg_infile->filename[infilenum],"r");
  
    if (infile == NULL) {
    
      fprintf(stderr,"writhecomp: Can't open %s.\n",arg_infile->filename[infilenum]);
      exit(1);
    
    }
  
    int err_num;
    char err_string[1024];
    size_t err_str_len = 1024;

    writhe_pline = plc_read(infile,&err_num,err_string,err_str_len);

    if (err_num != 0) {
    
      fprintf(stderr,"writhecomp: File %s does not appear to contain VECT data.\n",
	      arg_infile->filename[infilenum]);
      continue;
      
    }
  
    fclose(infile);
  
    /* signal the user that we're alive and have loaded a file ok */
  
    if (VERBOSE) {
    
      fprintf(stderr,"writhecomp 0.9\n");
      fprintf(stderr,"Loaded %d component, %d vertex polyline from %s.",
	      writhe_pline->nc,plc_num_verts(writhe_pline),arg_infile->filename[infilenum]);
      
    }

    if (arg_resolution->count != 0) {

      plCurve *swap;
      plc_spline *swap_spline;
      int *nv,j;
      bool ok;
      double *length;

      length = calloc(writhe_pline->nc,sizeof(double));
      nv = calloc(writhe_pline->nc,sizeof(int));
      plc_arclength(writhe_pline,length);

      for(j=0;j<writhe_pline->nc;j++) {nv[j] = ceil(arg_resolution->dval[0]*length[j]);}  
       
      swap_spline = plc_convert_to_spline(writhe_pline,&ok);
      swap = plc_convert_from_spline(swap_spline,nv);

      plc_free(writhe_pline);
      free(nv);
      free(length);
      writhe_pline = swap;

      printf(stderr,"Resampled to resolution %g.\n",arg_resolution->dval[0]);

    }
    
#ifdef PARI
    
  /* Check the stack size before we start to compute. */
  
  if((writhe_pline.verts - 2) * 12 > PARI_STACK_SIZE) {
    
    fprintf(stderr,"writhecomp: Short stack in PARI model.\n"\
	    "            writhe of %s cannot be computed.\n",
	    infile_name);
    fprintf(stderr,"            %s has %d verts - requiring %d "\
	    "bytes of stack vs. current %d...\n",
	    infile_name,
	    writhe_pline.verts,
	    (writhe_pline.verts - 2) * 12,
	    PARI_STACK_SIZE);
    fprintf(stderr,"writhecomp: Exiting.\n");
    exit(1);
    
  }
  
#endif
  
  // Check that we are performing the computation in enough precision.
  
  if (floor(log10((double)(plc_num_verts(writhe_pline)))) + 2 > DEC_PRECISION  ||
      DEC_PRECISION < 10) {
    
    printf("writhecomp: Low Precision! %s computation inaccurate.\n",
	   arg_infile->filename[infilenum]);
    printf("            %d digits are recommended for 2 place \n",
	   (int)(floor(log10((double)(plc_num_verts(writhe_pline))))+2));
    printf("            accuracy when computing writhe for a %d-vertex polyline.\n",
	   plc_num_verts(writhe_pline));
    printf("            Estimated accuracy (+/-) 10^(%d).\n",
	   (int)(DEC_PRECISION - floor(log10((double)(plc_num_verts(writhe_pline))))) );
    printf("\n");
    printf("            Use the --precision <n> flag to increase precision.\n");
  }
  
  // Parse the arithmetic model, and set up the appropriate hook functions.
  
  if (arithmodel == cdouble) {	/* C Double model */
    
    init_arith_model = &init_cdouble_model;
    kill_arith_model = &kill_cdouble_model;
    kill_np          = &free;
    
    double_to_np = &double_to_cdouble;
    np_to_double = &cdouble_to_double;
    np_to_string = &cdouble_to_string;
    
    np_add       = &cdouble_add;
    np_mul       = &cdouble_mul;
    np_div	 = &cdouble_div;
    
    np_acos      = &cdouble_acos;
    np_sqrt      = &cdouble_sqrt;
    np_leq       = &cdouble_leq;
    
  } else if (arithmodel == mapm) {
    
    init_arith_model = &init_mapm_model;
    kill_arith_model = &kill_mapm_model;
    kill_np          = &kill_mapm;
    
    double_to_np = &double_to_mapm;
    np_to_double = &mapm_to_double;
    np_to_string = &mapm_to_string;
    
    np_add       = &mapm_add;
    np_mul       = &mapm_mul;
    np_div	 = &mapm_div;
    
    np_acos      = &mapm_acos;
    np_sqrt      = &mapm_sqrt;
    np_leq       = &mapm_leq;
    
  } else if (arg_fast->count == 0) {
    
    fprintf(stderr,"writhecomp: This version not compiled with model %d. \n",
	    arithmodel);
    exit(1);
    
  }
  
  // We are now ready to compute.
  
  if (SAFE == TRUE && arg_fast->count == 0) {
    
    if (!arithmodel_ok()) {
      
      fprintf(stderr,"writhecomp: Failed arithmetic model sanity checks.\n");
      exit(1);
      
    }
    
  }

  if (arg_fast->count != 0) {

    pwrithe = fast_writhe(writhe_pline);

  } else if (arg_fast_spline->count != 0) {

    pwrithe = fast_spline_writhe(writhe_pline);

  } else {
  
    pwrithe = writhe(writhe_pline,DEC_PRECISION,NULL);

  }
  
  if (SAFE == TRUE) { // Compute again in higher and lower precision.
    
    lpwrithe = writhe(writhe_pline,DEC_PRECISION - 5,NULL);
    hpwrithe = writhe(writhe_pline,DEC_PRECISION + 5,NULL);
    
    if (fabs(pwrithe - hpwrithe) < 1e-8 && fabs(pwrithe - lpwrithe) < 1e-8) {
      
      /* If we pass the SAFE mode test, print a result and terminate normally. */
      
      printf("writhe: %16g \n\n",pwrithe);
      fprintf(stderr,"SAFE mode certifies these results as likely ok.\n");
      plc_free(writhe_pline);
       
    } else {
      
      /* If we fail the SAFE mode test, DON'T return a result to
	 stdout, since if we're running in batch mode, we WANT to crash
	 the calling script. DO return all three values to stderr, in
	 case a user is watching and ready to try again. Also,
	 terminate abnormally to give a calling script a chance to
	 recover gracefully. */
      
      fprintf(stderr,"writhe: %16g (lower precision: %16g) "\
	      "(higher precision: %16g) \n",pwrithe,lpwrithe,hpwrithe);
      fprintf(stderr,"SAFE mode cannot guarantee 8 place accuracy.\n"\
	      "Suggest you recompute with more precision. \n");
      plc_free(writhe_pline);
       
    }
    
  } else {
    
    /* In regular mode, just print the results, in script-friendly
       form (if we wanted more, we'd turn on VERBOSE mode, remember). */
    
    printf("writhe: %16g \n",pwrithe);
    plc_free(writhe_pline);
     
  }
  
  }

  exit(0);

}
