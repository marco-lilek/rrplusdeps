/********************************************************************
vect2kp.c : This program takes an argument

         vect2kp <file>

         and converts knot geometry in the VECT format
         to a format readable by Knotplot.  
	 We expect the input filename 
         to be in the form <name.vect>, and write the 
         output filename in the form <name.kp>.
  
  This code was modified from vect_to_fe.c by Matt Mastin,
  original code written by Ted Ashton.
********************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include<plCurve.h>

int main(int argc, char *argv[]){
  
  FILE *infile, *outfile;
  plCurve *L;
  int edges, components;
  double length;
  int error_num;
  char error_str[80];
  
  if (argc == 1){
    printf(
      "Usage: Converts knot geometry in the VECT format for use in knotplot.\n"
      "Sample command line: ./vect2kp ~/TightKnots/vects/[.vect file]\n");
    return 0;
  }
  
  if (argc > 2){
    printf("That is not a valid filename \n");
    return 0;
  }
  
  if (argc < 2){
    printf("vect2kp.c : This program takes a file name and an inradius\n\n"
           "vect2kp <file> <inradius>\n\n"
           "and converts the .vect file to a .kp, readable by Knotplot\n"
           "to a format readable by Knotplot. We expect the input filename\n"
           "to be in the form <name.vect>, and write the\n"
           "output filename in the form <name.kp>.\n");
  }
  
  else
    { 
      if (strstr(argv[1], ".vect") == NULL){
	printf("No VECT file specified.\n");
	exit(1);
      }
      
      /*Open the VECT file for reading*/
      
      infile = fopen(argv[1], "r");
      
      if (infile == NULL){
	printf("Could not open %s.\n", argv[1]);
	exit(1);
      }
      
      /*Generates an output filename from the input*/
      char filename[strlen(argv[1])-2];
      strncpy(filename, argv[1], strlen(argv[1])-2);
      
      char *ps = strstr(filename, ".ve");
      
      char *add = ".kp";
      
      while(*ps != '\0'){
	*ps++ = *add++;
      }
      
      /*Open the kp file for writing*/
      
      outfile = fopen(filename, "w");
      
      if(outfile == NULL){
	printf("Could not open %s.\n", argv[1]);
	exit(1);
      }
      
      /*Load a link*/
      
      L = plc_read(infile,&error_num,error_str,80);
      if (error_num != 0) {
	fprintf(stderr,"%s:%s",argv[0],error_str);
	exit(1);
      }
      edges = plc_num_edges(L);
      components = (*L).nc;
      length = plc_arclength(L,NULL);

      /* Add the vertices to the outfile*/
      
      int i = 1;
      int count = 0;
      while(count < components){

	fprintf(outfile, "%s %d %s %d %s\n", "Component ", count+1, " of ", components, ":");

	int j;
	for(j=0; j < (*L).cp[count].nv; j++){
	  fprintf(outfile, "%lf %lf %lf\n", (*L).cp[count].vt[j].c[0],
		  (*L).cp[count].vt[j].c[1], (*L).cp[count].vt[j].c[2]);
	  i++;
	}
	count++;
      }
      
      
      plc_free(L);
      fclose(infile);
      fclose(outfile);
      
      printf("File %s written.\n",filename);
      
    }
}

