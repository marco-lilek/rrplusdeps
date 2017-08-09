/**********************************************************
 * This program will convert a knot plot coords file into
 * a .vect file, compatable with Geomview
 *
 *
 * Matt Mastin
 * University of Georgia
 * Math Department
 * ********************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<plCurve.h>

int main(int argc, char *argv[]){

  plCurve *L;
  FILE *infile, *outfile;
  char buffer[80];
  float x, y, z;
  int nc=0, i=0, j=0;

  if (argc == 1) {

    printf("Usage: kp2vect file.kp\n");
    exit(0);

  }

  char filename[strlen(argv[1])+2];
  
  //We first open the input file
  infile = fopen(argv[1], "r");

  if(infile == NULL){
   printf("Invalid input file, .kp file expected as parameter.\n");
   return 0;
   }

  //Generate name for output file
  strncpy(filename, argv[1], strlen(argv[1])); 

  char *ps = strstr(filename, ".kp");
  char *add = ".vect";
  
  while(*add != '\0'){
    *ps++ = *add++;
    }

  *ps = 0;  //Be sure to terminate the string
  
  //Now open the outfile stream
  outfile = fopen(filename, "w");
  
  //First parse infile to count components
  while(fgets(buffer, 80, infile) != NULL){
    //count the components
    if(buffer[0]=='C'){
      nc++;      
      }
    }

  //Now count vertices in each component
  int nv[nc];

  rewind(infile);

  //Initialize the vertex counts  
  for(i=0;i<nc;i++)
    nv[i]=0;

  //We need this index to start "before" the first line
  i=-1;

  //Count the vertices
  while(fgets(buffer, 80, infile) != NULL){
    if(buffer[0]=='C')
      i++;
    else
      nv[i]++;
    }  

  rewind(infile);

  //Contruct plCurve
  bool open[nc];
  int cc[nc];
  int colors = 1;
  double temp[3];

  for(i=0;i<nc;i++){
    open[i]=false;
    cc[i]=1;
    }
  
  L = plc_new(nc, nv, open, cc);

  //Write vertex data to L
  fgets(buffer, 80, infile);
  for(i=0;i<nc;i++){
    fgets(buffer, 80, infile);
    for(j=0;j<nv[i];j++){
      sscanf(buffer, "%lf %lf %lf", &((*L).cp[i].vt[j].c[0]), &((*L).cp[i].vt[j].c[1]), &((*L).cp[i].vt[j].c[2]));
      fgets(buffer, 80, infile);
      }
   }
  
  plc_fix_wrap(L);

  //Print the .vect file
  plc_write(outfile, L);  

  printf("File %s written.\n", filename);
  plc_free(L);
  fclose(infile);
  fclose(outfile);
}
