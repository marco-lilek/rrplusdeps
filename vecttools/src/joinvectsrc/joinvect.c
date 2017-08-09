/*******************************************************************

 joinvect.c : Joins together VECT curve.
           See the man page for details.
	   
	   ************************************************************/

#include<joinvect.h>
#include<torusdistance.h>


/**********************************************************************/

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_lit  *fullauto;
struct arg_lit  *split;
struct arg_lit  *torus;
struct arg_lit  *greedy;

struct arg_end *end;
struct arg_end *helpend;

nplCurve     *core;
FILE         *infile_fptr,*outfile_fptr;
double        twidth = 1e10;

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,nerrors;
  char           *outfile_name;
 
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     fullauto = arg_lit0("a","automatic-mode","decides on join order"),
     split = arg_lit0("s","split","splits the components of joined curve into separate vect files"),
     greedy = arg_lit0("g","greedy","uses a greedy algorithm instead of the full travelling salesman to order curves"),
     torus = arg_lit0("t","torus_mode","computes the join in the flat torus T^n (expects width of torus in a file toruswidth.txt in cwd"),
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

      arg_print_errors(stdout,end,"joinvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("joinvect concatenates VECT files\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }
  
  /* Now we have parsed the arguments and are ready to work. */

  if (torus->count > 0) { 

    twidth = 1.0;

    FILE *infile;

    infile = fopen("toruswidth.txt","r");

    if (infile != NULL) {
      
      fscanf(infile,"%lf",&twidth);
      fclose(infile);
 
    } 

    printf("joinvect: Computing distances in torus of width %g.\n",twidth);

  }


  nplCurve **curves;
  curves = calloc(corefile->count,sizeof(nplCurve *));
  nplCurve *thiscurve;
 
  int ncurves = 0;
 
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"joinvect: Couldn't open file %s.\n",corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }
  
     int plr_error_num;
     char plr_error_str[1024];

     thiscurve = nplc_read(infile_fptr,
			   &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"joinvect: link reading error\n%s\n",plr_error_str);
       ncurves--;
       continue;  /* Try the next file */
    
     }
  
     if (thiscurve->nc > 1) {

       fprintf(stderr,
	       "joinvect: Curve %s has more than one component. \n"
	       "          Use splitvect to separate components \n"
	       "          before joining.\n",corefile->filename[infilenum]);

       nplc_free(thiscurve);
       continue;
     }

     if (thiscurve->cp[0].open == false) {

       fprintf(stderr,
	       "joinvect: Curve %s is already closed. Skipping it.\n",
	       corefile->filename[infilenum]);
       
       nplc_free(thiscurve);
       continue;
     }

     fclose(infile_fptr);

     /* Now add this guy to the list. */

     curves[ncurves++] = nplc_copy(thiscurve);
     nplc_free(thiscurve);

  }

  printf("joinvect: Loaded %d curves from VECT files.\n",ncurves);

  /* Now determine order. */

  int *order;
  order = calloc(ncurves,sizeof(int));
  int i;

  for(i=0;i<ncurves;i++) {

    order[i] = i+1;

  }

  order[ncurves-1]=0;

  if (fullauto->count > 0) { 

    printf("joinvect: Searching for optimal join order ");

    if (greedy->count > 0) {

      printf("using greedy algorithm.\n");
      greedy_order(ncurves,curves,order,twidth);

    } else {

      printf("using full TSP algorithm.\n");
      tsp_order(ncurves,curves,order,twidth);

    }

  }

  /* Now actually join the curves. */

  nplCurve *joined;
  bool    *open;
  int     *cc;
  bool    *used;

  /* Here's the rub... we may have more than one component here. */

  used = calloc(ncurves,sizeof(bool));

  int nc=0,startc,thisc;
  int *nv;

  /* First, determine the number of loops and the size of each. */

  nv = calloc(ncurves,sizeof(int));
  open = calloc(ncurves,sizeof(bool));
  cc = calloc(ncurves,sizeof(int));

  for(startc=0;startc<ncurves;startc++) {

    if (!used[startc]) {

      used[startc] = true;
      nv[nc] += nplc_num_verts(curves[startc]);
      
      for(thisc=order[startc];thisc != startc;thisc = order[thisc]) {
	
	nv[nc] += nplc_num_verts(curves[thisc]);
	used[thisc] = true;

      }

      nc++;

    }

  }

  printf("\njoinvect: Joined curves have %d components.\n",nc);
  joined = nplc_new(nplc_dim(curves[0]),nc,nv,open,cc);

  /* Now run the same loop, adding the actual vertices to the curves. */

  int cv;  
  nc = 0;
  int lastc;
  
  int long_cp,long_vt0,long_vt1,long_c;
  int cp,vt0,vt1;
  double long_len = -1,len;

  for(i=0;i<ncurves;i++) { used[i] = false; }

  for(startc=0;startc<ncurves;startc++) {
    
    if (!used[startc]) {

      printf("joinvect: Component loop %d-",startc);

      used[startc] = true;
      for(cv=0;cv<curves[startc]->cp[0].nv;cv++) {
	nplc_vect_copy(joined->cp[nc].vt[cv],curves[startc]->cp[0].vt[cv]);
      }

      /* Now search for longest edge. */
      long_c = startc;
      longest_edge(curves[startc],&long_cp,&long_vt0,&long_vt1,&long_len,twidth);
      
      for(lastc=startc,thisc=order[startc];thisc != startc;lastc=thisc,thisc = order[thisc]) {
	
	for(i=0;i<curves[thisc]->cp[0].nv;i++) {
	  nplc_vect_copy(joined->cp[nc].vt[cv],curves[thisc]->cp[0].vt[i]);
	  cv++;
	}
	used[thisc] = true;
	printf("(%3g)-%d-",torus_distance(nplc_end(curves[lastc],0),nplc_start(curves[thisc],0),twidth),thisc);

	/* Update longest edge calculation */
	longest_edge(curves[thisc],&cp,&vt0,&vt1,&len,twidth);
	if (len > long_len) {
	  long_c = thisc; long_cp = cp; long_vt0 = vt0; long_vt1 = vt1; long_len = len;
	}

      }
      
      printf("(%3g)-%d has %d vertices.\n",
	     torus_distance(nplc_end(curves[lastc],0),nplc_start(curves[startc],0),twidth),
	     startc,joined->cp[nc].nv);
      nc++;

      printf("joinvect: Longest edge %g at verts (%d-%d) of curve %d.\n\n",
	     long_len,long_vt0,long_vt1,long_c);
      
    }
    
  }

  printf("joinvect: Curve has average edgelength of %g.\n",nplc_arclength(joined,NULL)/nplc_edges(joined,NULL));
  printf("joinvect: If maximum join distance >> average edgelength, be careful!\n");
  
  nplc_fix_wrap(joined);

  free(open);
  free(nv);
  free(used);

  if (split->count == 0) { 
    
    /* Now output to a file */
    
    outfile_name = mangle(corefile->basename[0],".vect",".join.vect");
    outfile_fptr = fopen_or_die(outfile_name,"w");
    
    nplc_write(outfile_fptr,joined);
    fclose(outfile_fptr);
    nplc_free(joined);
    
    printf("joinvect: Wrote output file to %s.\n",outfile_name);
    free(outfile_name);

  } else {

    nplCurve *thiscp;
    int q;
    char cpname[1024];
    

    for(q=0;q<joined->nc;q++) {

      thiscp = nplc_component(q,joined);
      sprintf(cpname,"joined.%03d.vect",q);
      outfile_fptr = fopen_or_die(cpname,"w");
      
      nplc_write(outfile_fptr,thiscp);
      fclose(outfile_fptr);
      printf("joinvect: Wrote %d vertex loop to output file to %s.\n",
	     thiscp->cp[0].nv,cpname);

      nplc_free(thiscp);
    
    }

    nplc_free(joined);

  }
  
  /* Now do housekeeping */

  for(i=0;i<ncurves;i++) {
    
    nplc_free(curves[i]);
    
  }
  
  exit(0);
  
}




