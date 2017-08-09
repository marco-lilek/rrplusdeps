/*******************************************************************

   whittenvect.c : Applies Whitten group operations to a plCurve
                   given as a vect file. 
                   See the man page for details.
	   
	   ************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef STRING_H
#include <string.h>
#endif

#include<plCurve.h>
#include<../utilib/mangle.h>
#include<../utilib/ordie.h>
#include<argtable2.h>

/* Global variables live here. */

struct arg_str  *whittenop;
struct arg_lit  *verbose;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_lit  *quiet;
struct arg_file *arg_outfile;

struct arg_end *end;
struct arg_end *helpend;

plCurve *core;
FILE    *infile_fptr,*outfile_fptr;

int    QUIET=0;

/****************************** Main procedure ********************************/

/* Calls strtok and fails on error */
#define mystrtok(A,B)	\
  thistok = strtok((A),(B)); \
  if (thistok == NULL) {					\
     printf("whittenvect: Whitten operation %s should have %d " \
	    "epsilon elements", whittenop, nc); \
     exit(1); } \
  thistok = thistok

#define myeps(A,B) \
   if (sscanf((A),"%d",&(B)) != 1) { \
     printf("whittenvect: Couldn't parse element %s in "\
	    "whitten operation %s as integer.\n", (A),whittenop); \
     exit(1); } \
   if ((B) != 1 && (B) != -1) {					\
     printf("whittenvect: Couldn't parse element %s in "\
	    "whitten operation %s as +1 or -1.\n", (A), whittenop); \
     exit(1); } \
   thistok = thistok			   

void parse_whittenop(int nc,const char *whittenop,int *mirror,int **epsilon,int **perm)
/* Parse a whitten operation in standard form into our special internal form */
{ 

  char *opcopy;

  opcopy = calloc(strlen(whittenop)+10,sizeof(char));
  strcpy(opcopy,whittenop);

  /* Eliminate leading and trailing parentheses. */

  if (opcopy[0] != '(') {

    printf("whittenvect: Whitten operation %s must start with (\n",whittenop);
    exit(1);

  } else { 

    opcopy[0] = ' ';

  }

  int oplen;
  oplen = strlen(opcopy);
  
  if (opcopy[oplen-1] != ')') {

    printf("whittenvect: Whitten operation %s must end with )\n",whittenop);
    exit(1);

  } else { 

    opcopy[oplen-1] = 0;

  }

  /* Tokenize on ',' to separate e_i and permutation */

  int  i;
  char *thistok;

  *perm = calloc(2*nc,sizeof(int)); /* This is too big, but that's ok. */
  *epsilon = calloc(nc,sizeof(int)); 

  mystrtok(opcopy,",");
  myeps(thistok,(*mirror));

  for(i=0;i<nc;i++) {

    mystrtok(NULL,",");  // Get subsequent tokens from opcopy string
    myeps(thistok,((*epsilon)[i]));

  }
  
  /* The next token should be the permutation, which is either a bunch of 
     cycles or the character 'e'. Our internal format for permutations is 
     an array of 2n integers so that the integer after the first occurence of
     a given integer is its image under the permutation. Integers which 
     do not occur in the string are not changed. 

     We also change from 1 based notation to 0 based notation.

     Example: 3 elements (1 2) -> 0 1 0 -1 -1 -1. 
     Example: 5 elements (1 3 2)(4 5) -> 0 2 1 0 3 4 3 -1 -1 -1. */

  int permn = 0;

  mystrtok(NULL,",");  // Again, we're reading from opcopy
  if (strchr(thistok,'e') != NULL) {  

    for(i=0;i<2*nc;i++) {
      (*perm)[i] = -1; /* No elements occur in the permutation */
    }

  } else {
    
    /* We now rewrite the permutation in the form above */
    /* Make a copy to re-tokenize on () to split into cycles. */

    char *permcopy,*thiscycle;
    char **cycles;
    int  ncycles = 0,i;

    permcopy = calloc(strlen(thistok)+10,sizeof(char));
    strcpy(permcopy,thistok);
    cycles = calloc(nc,sizeof(char *));
    
    for(thiscycle = strtok(permcopy,"()");
	thiscycle != NULL;
	thiscycle = strtok(NULL,"()")) {

      /* Copy the cycle into the cycles array. */

      cycles[ncycles] = calloc(strlen(thiscycle)+10,sizeof(char));
      strcpy(cycles[ncycles],thiscycle);
      ncycles++;

    }

    if (verbose->count > 0) {

      printf("Found %d cycles in whitten element %s.\n",ncycles,whittenop);
      
    }

    /* Now retokenize each cycle */

    for (i=0;i<ncycles;i++) {

      char *thiselt;
      int   cstart = permn; /* Position of first elt in cycle */

      /* Now tokenize on " " to split into individual indices */
      
      for( thiselt = strtok(cycles[i]," ");
	   thiselt != NULL;
	   thiselt = strtok(NULL," ") ) {
	
	if (sscanf(thiselt,"%d",&((*perm)[permn])) != 1) {
	  
	  printf("whittenvect: Could not parse cycle element %s in %s.\n",
		 thiselt,whittenop);
	  exit(1);
	  
	}
	
	(*perm)[permn]--;  /* Switch to 0-based indexing */
	permn++;

      }
	 
      /* We have reached the end of the cycle, so recopy the first elt. */
      /* This is already 0-based, so we don't need to fix it here. */

      (*perm)[permn++] = (*perm)[cstart];
      free(cycles[i]); /* We don't need cycles[i] anymore */
      
    } /* end of cycles */

    free(cycles); /* Don't need this guy anymore */

    /* We now fill the rest of the perm array with -1s in order to be safe */

    for(;permn < 2*nc;permn++) { (*perm)[permn] = -1; }
    free(permcopy);

  } /* end of cycle parsing */

  /* We should now have tokenized completely. */

  if (verbose->count > 0) { /* Report to the user. */
    
    int i;
    
    printf("Converted whitten operation %s into:\n",whittenop);
    printf("  Mirror:   %d.\n",*mirror);

    printf("  Epsilons: ");
    for(i=0;i<nc;i++) { printf("%d ",(*epsilon)[i]); }
    printf("\n");

    printf("  Permutation: ");
    for(i=0;i<2*nc;i++) { printf("%d ",(*perm)[i]); }
    printf("\n");

  }

  /* Now free the various stupid strings we've allocated. */

  free(opcopy);

}

char *maketail(const char *whittenop)
/* Converts the whittenop to a filename-friendly string by replacing */
/* whitespace with -. */
{
  char *tail;
  char *ptr;

  tail = calloc(strlen(whittenop)+40,sizeof(char));
  sprintf(tail,".%s.vect",whittenop);
  
  for(ptr=tail;*ptr!=0;ptr++) { 

    if (*ptr == ' ') { *ptr = '-'; } 
    if (*ptr == '(') { *ptr = '['; }
    if (*ptr == ')') { *ptr = ']'; }

  };


  return tail;
}
  
  
int main(int argc,char *argv[])
{
  int         infilenum,nerrors;
  
  void *argtable[] = 
    {
     whittenop = arg_strn("o","operation","(1,-1,-1,(2 1))",
			  1,100000,"group operation"),
     verbose = arg_lit0("v","verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     arg_outfile = arg_filen(NULL,"output","<file>",0,100000,"output files"),
     quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};

  fprintf(stdout,"whittenvect (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("whittenvect: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"whittenvect");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("whittenvect applies whitten group operation to VECT file.\n"
	     "\n"
	     "The Whitten group acts on a link of n components by mirroring,\n"
	     "reversing orientations of components and permuting components.\n"
	     "An element of the group is written in the form:\n"
	     "\n"
	     "(+-1, +-1, ... , +-1, p)\n"
	     "\n"
	     "where the first +-1 decides whether the link is mirrored, and\n"
	     "the next n +-1s determine the orientations of the components \n"
	     "according to the rule\n"
	     "\n"
	     "(e_0,e_1,...e_n,p)(K_1,...,K_n) = (e_1 K_p(1), ..., e_n K_p(n))\n"
	     "\n"
	     "The permutation p must be written as 'e' (for the identity) or\n"
	     "in cycle form with the cycles separated by parentheses\n"
	     "\n"
	     "(1 3 5)(2 4), or (1 3)(7 13)(2 3)\n"
	     "\n"
	     "the elements of the cycles must be disjoint.\n\n"
	     "\n"
	     "Example: whittenvect -o (1,-1,1,(1 2)) hopflink.vect\n"
	     "         whittenvect -o (-1,1,e) 3_1.vect\n"
	     "\n"
             "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  QUIET = quiet->count > 0; /* Register if we're in batch mode */
  
  /* Now we have parsed the arguments and are ready to work. */
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We begin by loading the current core link. */

     infile_fptr = fopen(corefile->filename[infilenum],"r");
  
     if (infile_fptr == NULL) {
    
       fprintf(stderr,"whittenvect: Couldn't open file %s.\n",
	       corefile->filename[infilenum]);
       continue;  /* Try the next file */
 
     }

     int plr_error_num;
     char plr_error_str[1024];
  
     core = plc_read(infile_fptr,
		     &plr_error_num,plr_error_str,sizeof(plr_error_str));
  
     /* We now demonstrate the plCurve library's error handling protocol: */
  
     if (plr_error_num > 0) {   /* This is the signal for an error. */
    
       fprintf(stderr,"whittenvect: link reading error\n%s\n",plr_error_str);
       continue;  /* Try the next file */
    
     }
  
     fclose(infile_fptr);

     /* Report to the user */

     if (!QUIET) {

       printf("%s\n",corefile->basename[infilenum]);

     }

     /* Now we actually do the whitten op. */

     int mir,*eps,*perm;

     parse_whittenop(core->nc,whittenop->sval[infilenum % whittenop->count],
		     &mir,&eps,&perm);

     plc_whitten(core,mir,eps,perm);

     /* Now we rename */

     char *newtail,*ofname;
     FILE *outfile;

     if (arg_outfile->count > 0) {

       outfile = fopen(arg_outfile->filename[infilenum%arg_outfile->count],"w");
       if (outfile == NULL) {

	 fprintf(stderr,"whittenvect: Couldn't open output file %s\n",arg_outfile->filename[infilenum%arg_outfile->count]);
	 exit(1);

       }

       plc_write(outfile,core);

     } else {

       newtail = maketail(whittenop->sval[infilenum % whittenop->count]);
       outfile = fmangle(corefile->filename[infilenum],".vect",newtail);
       plc_write(outfile,core);
       
       if (!QUIET) {
       
	 ofname = mangle(corefile->filename[infilenum],".vect",newtail);
	 printf("Wrote modified VECT file to: %s\n",ofname);
	 free(ofname);

       }

     }
     
     /* Now we cleanup memory. */

     plc_free(core);

  }

  return 0;

}




