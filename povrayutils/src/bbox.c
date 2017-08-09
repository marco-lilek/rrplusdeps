#include "bbox.h"

void write_bbox(plc_vector l,plc_vector u,FILE *outfile)
{
  fprintf(outfile,"BBOX\n");
  fprintf(outfile,"%g %g %g\n%g %g %g\n",
	  plc_M_clist(l),plc_M_clist(u));
}
  
bool read_bbox(plc_vector *l,plc_vector *u,FILE *infile)
{
  if (fscanf(infile,"BBOX %lf %lf %lf %lf %lf %lf ",
	     &(l->c[0]),&(l->c[1]),&(l->c[2]),
	     &(u->c[0]),&(u->c[1]),&(u->c[2])) != 6) return false;

  return true;
}
	 


  
