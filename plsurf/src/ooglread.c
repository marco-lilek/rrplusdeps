/************************************************/
/*             Oogl reading procedures          */
/************************************************/

/* This code is part of the linked list library.*/

#include "plsurf.h"
#include "plsurf_internal.h"

int skip_whitespace_and_comments(FILE *infile)

     /* Procedure positions the file pointer on */
     /* next non-whitespace character, returning */
     /* FALSE if EOF happens first. We skip anything */
     /* between a # and a newline. */

{
  int thischar,commentflag = {0};

  /* First, we check to make sure that infile looks legit. */

  if (infile == NULL) {

    fprintf(stderr,"skip_whitespace_and_comments: infile is a null pointer.\n");
    exit(2);
    
  }
  
  /* Now we start to work. */

  for(;;) {

    thischar = fgetc(infile);

    if (thischar == EOF) {  /* Reached end of file before a non-space, non-comment */
      
      return 0;
    
    } else if (thischar == '#') { /* Started a comment. */
    
      commentflag = 1;

    } else if (thischar == '\n' && commentflag) { /* End a comment. */

      commentflag = 0;

    } else if (!isspace(thischar) && !commentflag) { /* Found a hit! */

      ungetc(thischar,infile);
      return 1;

    } /* It must have been a space or a non-space in a comment. */

  }

}

int scandoubles(FILE *infile,int ndoubles, ...)

     /* Procedure scans for nfloats floating point (or double) numbers, ignoring 
	whitespace and comments between them. We expect the variable length arguments 
	to contain a collection of pointers to doubles. If not, there's trouble. */

{

  int nconverted = 0,i;
  va_list ap;
  double *thisdouble;

  /* First, we check for overall sanity. */

  if (infile == NULL) {

    fprintf(stderr,"scandoubles: infile is a null pointer.\n");
    exit(2);

  }

  if (ndoubles < 1) {

    fprintf(stderr,"scandoubles: ndoubles (%d) is less than one.\n",ndoubles);
    exit(2);

  }

  va_start(ap,ndoubles);

  /* Now we're ready to work. */

  for (i=0;i<ndoubles;i++) {	/* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */

      return nconverted;

    }

    thisdouble = va_arg(ap,double *);

    if (fscanf(infile,"%lf",thisdouble) != 1) {	/* We couldn't scan. */

      return nconverted;	/* So give up here */

    } else {			/* Else record our victory */

      nconverted++;

    }

  }

  va_end(ap);

  return nconverted;

}

int scanints(FILE *infile,int nints, ...)

     /* Procedure scans for nints integers, ignoring whitespace and
	comments between them. We expect the variable length arguments
	to contain a collection of pointers to doubles. If not,
	there's trouble. */

{

  int nconverted = 0,i;
  va_list ap;
  int *thisint;

  /* First, we check for overall sanity. */

  if (infile == NULL) {

    fprintf(stderr,"scanints: infile is a null pointer.\n");
    exit(2);

  }

  if (nints < 1) {

    fprintf(stderr,"scanints: nints (%d) is less than one.\n",nints);
    exit(2);

  }

  va_start(ap,nints);

  /* Now we're ready to work. */

  for (i=0;i<nints;i++) {	/* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */

      return nconverted;

    }

    thisint = va_arg(ap,int *);

    if (fscanf(infile,"%d",thisint) != 1) {	/* We couldn't scan. */

      return nconverted;	/* So give up here */

    } else {			/* Else record our victory */

      nconverted++;

    }

  }

  va_end(ap);

  return nconverted;

}

int skipwhitespace(FILE *infile)

     /* Procedure skips whitespace (but NOT newlines). */
     /* Returns 1 if we find something, 0 if EOF. */
{
  int thischar;

  for(;;) {

    thischar = fgetc(infile);
    
    if (thischar == '\n' || !isspace(thischar)) {

      ungetc(thischar,infile);
      return 1;

    } else if (thischar == EOF) {

      return 0;

    }

  }

}

int scancolor(FILE *infile,plc_color *thiscolor)

     /* Procedure reads a color in Geomview RGBA format. */
     /* The only tricky part is that the color could be  */
     /* RGB (doubles) or RGB (ints between 0 and 255) or */
     /* RGBA in either format. */

     /* The trick here is that color information is line-based */
     /* so we can't use our ordinary scanning stuff. */

     /* We return 0 if we weren't able to scan a color, */
     /* 1 otherwise. The file will be at the newline or # */
     /* marking the end of this line. */

{
  int thischar;

  /* First, we make sure that the pointers look ok. */

  if (infile == NULL || thiscolor == NULL) {

    fprintf(stderr,"scancolor: Infile or thiscolor is a null pointer.\n");
    exit(2);

  }

  /* Now we go ahead and work. */

  if (!skipwhitespace(infile)) return 0; /* Skip ahead. */

  thischar = fgetc(infile);	         /* Test for '\n' */
  if (thischar == '\n' || thischar == '#') return 0;       
  ungetc(thischar,infile);

  if (fscanf(infile,"%lf",&(thiscolor->r)) != 1) return 0; /* Read the number. */

  /* Now read G */

  if (!skipwhitespace(infile)) return 0; /* Skip ahead. */

  thischar = fgetc(infile);	         /* Test for '\n' */
  if (thischar == '\n' || thischar == '#') return 0;       
  ungetc(thischar,infile);

  if (fscanf(infile,"%lf",&(thiscolor->g)) != 1) return 0; /* Read the number. */

  /* Now read B */

  if (!skipwhitespace(infile)) return 0; /* Skip ahead. */

  thischar = fgetc(infile);	         /* Test for '\n' */
  if (thischar == '\n' || thischar == '#') return 0;       
  ungetc(thischar,infile);

  if (fscanf(infile,"%lf",&(thiscolor->b)) != 1) return 0; /* Read the number. */

  /* A is special, as it is allowed to not be present. */

  if (!skipwhitespace(infile)) return 0; /* Skip ahead. */

  thischar = fgetc(infile);	         /* Test for '\n' */
  
  if (thischar != '\n' && thischar != '#') { /* There is an "A" */      

    ungetc(thischar,infile);
    if (fscanf(infile,"%lf",&(thiscolor->alpha)) != 1) return 0; /* Read the number. */

  } else {

    thiscolor->alpha = 1.0;

  }

  /* We have now read all 3 or 4 numbers. We next check whether they */
  /* are integers or floating point values. If they are integers, they */
  /* should be in the range 0..255. If they are 0 or 1, the color */
  /* is ambiguous, and we choose to interpret it in the floating point sense. */

  if (thiscolor->r > 1.5) thiscolor->r /= 255.0;
  if (thiscolor->g > 1.5) thiscolor->g /= 255.0;
  if (thiscolor->b > 1.5) thiscolor->b /= 255.0;
  if (thiscolor->alpha > 1.5) thiscolor->alpha /= 255.0;
  
  /* Now we check for sanity. */

  if (!goodcolor(*thiscolor)) return 0;

  return 1;

}
