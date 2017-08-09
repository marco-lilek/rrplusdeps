/*******************************************************************

mkellipse.c
	   
	   ********************************************************/

#include<orientvect.h>
#include<../utilib/mangle.h>
#include<../utilib/ordie.h>


int main()
{
  int nv = {100};
  int cc = {0};
  bool open = {false};

  plCurve *L;

  L = plc_new(1,&nv,&open,&cc);

  double theta;
  int i;

  for(i=0,theta=0;i<100;i++,theta+=6.28318531/100) {
    
    L->cp[0].vt[i] = plc_build_vect(2*cos(theta),sin(theta),0);

  }

  FILE *ell;

  ell = fopen("ellipse.vect","w");
  plc_write(ell,L);
  fclose(ell);
  plc_free(L);

  printf("Wrote ellipse to ellipse.vect\n");

  return 0;

}
