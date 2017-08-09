#include "tube.h"
#include "uv.h"

void write_uvfile(struct uvbuf *uvb,FILE *outfile)

     /* Writes the uv buffer to a simple text file. */
     /* The format is 

     UVDATA
     nverts
     u_0 v_0
     u_1 v_1
     ...
     u_nverts-1 v_nverts-1

     */

{
  assert(outfile != NULL);
  assert(uvb != NULL);

  fprintf(outfile,
	  "UVDATA\n"
	  "%d\n",uvb->verts);

  int i;
  
  for(i=0;i<uvb->verts;i++) {

    fprintf(outfile,"%g %g\n",uvb->uv[i].c[0],uvb->uv[i].c[1]);
    
  }

}
  

int  load_uvfile(struct uvbuf *uvb,FILE *infile)

     /* Reads from a simple format text file. Will allocate 
        storage in uvbuf as needed. Returns 1 on success, 0 on error. */

{
  assert(infile != NULL);
  assert(uvb != NULL);

  if (fscanf(infile,
	     "UVDATA\n"
	     "%d\n",&(uvb->verts)) != 1) {

    return 0;

    }

  assert(uvb->verts > 0);

  uvb->uv = calloc(uvb->verts,sizeof(plc_vector));

  int i;
  
  for(i=0;i<uvb->verts;i++) {

    if (fscanf(infile,"%lf %lf\n",&(uvb->uv[i].c[0]),&(uvb->uv[i].c[1])) != 2) {

      return 0;
    
    }
  }  

  return 1;

}
