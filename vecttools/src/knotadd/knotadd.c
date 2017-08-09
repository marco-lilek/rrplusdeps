#include<config.h>

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_MATH_H
#include <math.h>
#endif

//#include <octrope.h>

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#include <plCurve.h>
#include <argtable2.h>

//If the biggest splitting radius is smaller than this
//then we've screwed the moose.
//
//have you ever seen a moose??? they're quite large... i don't know if they'd notice...

double THRESHOLD = .00001;

//copies one vector to another

void vect_copy (plc_vector * v1, plc_vector * v2) {
  
  v1->c[0] = v2->c[0];
  v1->c[1] = v2->c[1]; 
  v1->c[2] = v2->c[2];

  return;

}

//calculates the distance between projections of two
//vectors in the xy-plane.

double twod_len(plc_vector *Point1, plc_vector *Point2)
{
  plc_vector Vector;

  Vector.c[0] = Point2->c[0] - Point1->c[0];
  Vector.c[1] = Point2->c[1] - Point1->c[1];

  return (double)sqrt(Vector.c[0]*Vector.c[0] + Vector.c[1]*Vector.c[1]);
}

//calculates the distance between two points.

double threed_len(plc_vector *Point1, plc_vector *Point2)
{
  plc_vector Vector;

  Vector.c[0] = Point2->c[0] - Point1->c[0];
  Vector.c[1] = Point2->c[1] - Point1->c[1];
  Vector.c[2] = Point2->c[2] - Point1->c[2];
  
  return (double)sqrt(Vector.c[0]*Vector.c[0] + Vector.c[1]*Vector.c[1] + Vector.c[2]*Vector.c[2]);
}

//calculates the orthogonal distance between point and the vector v2-v1.

double min_dist(plc_vector * v1, plc_vector * v2, plc_vector * point) {
  double  LineMag;
  double  U;
  plc_vector Intersection;
  double len1, len2;


  LineMag = twod_len(v1, v2);
  len1 = twod_len(point, v1);
  len2 = twod_len(point, v2);
 

  U = ( ( ( point->c[0] - v1->c[0] ) * ( v2->c[0] - v1->c[0] ) ) +
        ( ( point->c[1] - v1->c[1] ) * ( v2->c[1] - v1->c[1] ) ) ) /
    ( LineMag * LineMag );
  
  if( U < 0.0 || U > 1.0 ) {
    if(len1 < len2) {
      return len1;
    } else {
      return len2;
    }
  }

  Intersection.c[0] = v1->c[0] + U * ( v2->c[0] - v1->c[0] );
  Intersection.c[1] = v1->c[1] + U * ( v2->c[1] - v1->c[1] );
  Intersection.c[2] = 0;

  return twod_len( point, &Intersection );
}

//calculates the highest and lowest points on a link.

plc_vector * high_low (plCurve * link) {

  plc_vector * tb_points = calloc(2,sizeof(plc_vector));

  int x,y;

  for (x = 0; x < link->nc;x++) {
    for (y = 0; y < link->cp[x].nv; y++) {
      if (x == 0 && y == 0) {
	vect_copy(&tb_points[0],&(link->cp[x].vt[y]));
	vect_copy(&tb_points[1],&(link->cp[x].vt[y]));
      } else {
	if(link->cp[x].vt[y].c[2] < tb_points[0].c[2]) {
	  vect_copy(&tb_points[0],&(link->cp[x].vt[y]));
	}
	if(link->cp[x].vt[y].c[2] > tb_points[1].c[2]) {
	  vect_copy(&tb_points[1],&(link->cp[x].vt[y]));
	}
      }
    }
  }
  return tb_points;
}

struct arg_file *input_files;
struct arg_file *output_file;
struct arg_int  *components;
struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;
struct arg_end  *end, *helpend;

int main (int argc, char *argv[]) {

  int  nerrors;
  
  void *argtable[] = 
    {
      input_files = arg_filen(NULL,NULL,"<file>",2,2,"input files"),
      output_file = arg_file0("o","output","<file","output file"),
      components = arg_intn("c","component","<int>",0,2,"components to attach"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("knotadd: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"knotadd (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"knotadd");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"knotadd (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      printf("knotadd joins knots or links to make composites.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  //  printf("here");
  plCurve ** link;
  plCurve * link_join;

  

  FILE * infile1 = fopen(input_files->basename[0],"r");
  FILE * infile2 = fopen(input_files->basename[1],"r");
  FILE * outfile;

  if (output_file->count > 0) {

    outfile = fopen(output_file->basename[0],"w");

  } else {

    outfile = fopen("composite.vect","w");

  }

  if (infile1 == NULL || infile2 == NULL || outfile == NULL) {

    fprintf(stderr,"knotadd: Couldn't open input or output files.\n");
    exit(1);
    
  }

  int cn[2];  
  int w,x,y,z;

  int new_comp;
  int * new_nv;
  int * new_acyclic;
  int * new_colors;

  plc_vector * high_low1;
  plc_vector * high_low2;

  double trans_dist[3];
  double min_dist_temp;
  double min_dist_temp2;
  double height;

  double sel1_dist = 0;
  double sel2_dist = 0;

  int min_dist_vert[2];
  double min_dist_fin[2];
  double avg_dist[2];
  double average;
  int new_segs = 0;

  int flag; 
  int dflag;

  //argc[4] and argc[5], the components to add

  if (components->count > 0) {

    cn[0] = components->ival[0];

    if (components->count > 1) {

      cn[1] = components->ival[1];

    } else {

      cn[1] = 0;

    }

  } else {

    cn[0] = cn[1] = 0;

  }

  //read in the files
  
  int e0 = 0;
  int e1 = 0;
  
  int *e0_p = &e0;
  int *e1_p = &e1;
  
  char c1[5];
  char c2[5];
  
  size_t a=5;
  size_t b=5;
	
  link = calloc(2, sizeof(plCurve *));
  link[0] = plc_read(infile1,e0_p,c1,a);
  link[1] = plc_read(infile2,e1_p,c2,b);

  //check error number from plc_read and exit if error occurs

  if(e0!=0 || e1!=0){
    printf("%i error reading file",e0);
    return(0);  
  }

  for (w = 0; w < 2; w++) {
    avg_dist[w]=0;
    for(x = 0; x < link[w]->cp[cn[w]].nv; x++) {

      //here we start calculating the average segment length
      //on each link to see how to divide up the bridging segments
      //into shorter average segments

      avg_dist[w] += threed_len(&(link[w]->cp[cn[w]].vt[x]),&(link[w]->cp[cn[w]].vt[x+1]));

      //Label 1: here we loop through the vect in question, and for every vertex on the
      //component we're going to join (cn[w]), we see what the minimum distance is 
      //between it and other segments on the entire link. Whatever vertex has the 
      //biggest minimum distance is stored in min_dist_vert and min_dist_fin.
      //This vertex will give us the most room for splicing in the bridge (see
      //documentation).

      flag = 0;
      for(y = 0; y < link[w]->nc; y++) {
	for (z = 0; z < link[w]->cp[y].nv; z++) {
	  if((cn[w] == y) &&
	     (z == x || z+1 == x || (z+1 == link[w]->cp[cn[w]].nv && x == 0))) { 
	    continue;
	  }
	  if(flag == 0) {
	    min_dist_temp = min_dist(&(link[w]->cp[y].vt[z]),&(link[w]->cp[y].vt[z+1]),&(link[w]->cp[cn[w]].vt[x]));
	    flag = 1;
	  } else {
	    min_dist_temp2 = min_dist(&(link[w]->cp[y].vt[z]),&(link[w]->cp[y].vt[z+1]),&(link[w]->cp[cn[w]].vt[x]));
	    if (min_dist_temp2 < min_dist_temp) {
	      min_dist_temp = min_dist_temp2;
	    }
	  }
	}
      }
      if (x == 0 || min_dist_temp > min_dist_fin[w]) {
	min_dist_vert[w] = x;
	min_dist_fin[w] = min_dist_temp;
      }
    }

    avg_dist[w] = avg_dist[w]/link[w]->cp[cn[w]].nv;
  }
  //now we finish our calculation of average segment length
  average = (avg_dist[0]+avg_dist[1])/2;

  //whichever min_dist_fin is smaller,
  //we're bound to this circle for splitting
  // our vertices apart.

  if (min_dist_fin[0] < min_dist_fin[1]) {
    min_dist_fin[1] = min_dist_fin[0];
  } else {
    min_dist_fin[0] = min_dist_fin[1];
  }

  //If out min_dist_fin is less than threshold,
  //whatever that's set at, then the distance
  // is too small to do anything with without
  // evolver barfing.

  if (min_dist_fin[0] < THRESHOLD) {
    printf("Can't find wide enough radius for joining\n");
    exit(1);
  }

  //calculate the high and low points on each vect.

  high_low1 =  high_low(link[0]);
  high_low2 =  high_low(link[1]);

  //Label 2: In this next bit of code, we calculate what
  //we have to add in the x, y, and z direction to put
  //one vertex that we are going to split directly over
  //the other one in the z direction without any risk of
  //the links bumping uglies.

  if (high_low1[0].c[2] <= high_low2[1].c[2]) {

    trans_dist[2] = high_low2[1].c[2] - high_low1[0].c[2] + 1;

  } else {

    trans_dist[2] = 1;

  }
  
  if(link[0]->cp[cn[0]].vt[min_dist_vert[0]].c[0] < link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[0]) {
    trans_dist[0] = (double)fabs(link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[0] - link[0]->cp[cn[0]].vt[min_dist_vert[0]].c[0]);
  } else {
    trans_dist[0] = -1 * (double)fabs(link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[0] - link[0]->cp[cn[0]].vt[min_dist_vert[0]].c[0]);
  }
  
  if(link[0]->cp[cn[0]].vt[min_dist_vert[0]].c[1] < link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[1]) {
    trans_dist[1] = (double)fabs(link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[1] - link[0]->cp[cn[0]].vt[min_dist_vert[0]].c[1]);
  } else {
    trans_dist[1] = -1 * (double)fabs(link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[1] - link[0]->cp[cn[0]].vt[min_dist_vert[0]].c[1]);
  }
  
  //link[0] now gets translated over link[1].

  for (x = 0; x < link[0]->nc; x++) {
    for (y = 0; y < link[0]->cp[x].nv; y++) {
      link[0]->cp[x].vt[y].c[0] += trans_dist[0];
      link[0]->cp[x].vt[y].c[1] += trans_dist[1];
      link[0]->cp[x].vt[y].c[2] += trans_dist[2];
    }
  }

  //What's the distance between the two vertices we're splitting?
  height = fabs(link[0]->cp[cn[0]].vt[min_dist_vert[0]].c[2] - link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[2]);
  new_segs= ceil(height/average);
  average = height/new_segs;
  //we need new_segs worth of average length segments in order to bridge out knots
  //appropriately.


  //Label 3: We create a shell of a link to insert
  //our joined knot into.

  new_comp = link[0]->nc + link[1]->nc -1;
  new_nv = calloc(new_comp, sizeof(int));
  new_acyclic = calloc(new_comp, sizeof(int));
  new_colors = calloc(new_comp, sizeof(int));
  
  for (x = 0; x < new_comp; x++) {
    new_acyclic[x] = 0;
    new_colors[x] = 1;
  }
  //here we figure out just how many 
  //vertices each component of our joined knot will have.

  for (x = 0; x < link[0]->nc; x++) {
    new_nv[x] = link[0]->cp[x].nv;
  }
  for (x = 0; x < link[1]->nc; x++) {
    if(x == cn[1]) {
      if(argc == 7) {
	if(strcmp(argv[6],"-d") == 0){
	  new_nv[cn[0]] += link[1]->cp[x].nv + (new_segs*2) - 2; //To account for all the new vertices we're adding.
	  printf("Doubled vertices are silly, but I'll see what I can do.\n");
	  dflag = 1;
	} else {
	  printf("Final command line argument -- %s -- not understood.\n", argv[6]);	  
	  exit(1);
	}
      } else {
	new_nv[cn[0]] += link[1]->cp[x].nv + (new_segs*2); //To account for all the new vertices we're adding.
	dflag = 0;
      }
    } else {
      if(x < cn[1]) {
	new_nv[x+link[0]->nc] = link[1]->cp[x].nv;
      } else {
	new_nv[x + link[0]->nc - 1] = link[1]->cp[x].nv;
      }
    }
  }

  //Voila! It's created.
  link_join = plc_new(new_comp,new_nv,(bool*) new_acyclic,new_colors);
  
  //Label 4: Here's the big part. We copy in all the vertex info, both
  //unaffected and unaffected components.

  for (x = 0; x < link[0]->nc; x++) {//copying in unaffected components
    if (x == cn[0]) {
      continue;
    } 
    for (y = 0; y < link[0]->cp[x].nv; y++) {
      vect_copy(&(link_join->cp[x].vt[y]),&(link[0]->cp[x].vt[y]));
    }
  }

  for (x = 0; x < link[1]->nc; x++) {//copying in unaffected components
    if (x == cn[1]) {
      continue;
    } 
    for (y = 0; y < link[1]->cp[x].nv; y++) {
      if(x < cn[1]) {
	vect_copy(&(link_join->cp[x+link[0]->nc].vt[y]),&(link[1]->cp[x].vt[y]));
      } else {
	vect_copy(&(link_join->cp[x + link[0]->nc - 1].vt[y]),&(link[1]->cp[x].vt[y]));
      }
    }
  }

  //copying joined components...Here we go:

  //We start by creating one of the two new split vertices on
  //link[0]. We do this by copying the original guy and then moving him
  //rad/sqrt(2) over in the x = y direction (see documentation).

  vect_copy (&(link_join->cp[cn[0]].vt[0]),&(link[0]->cp[cn[0]].vt[min_dist_vert[0]]));
  link_join->cp[cn[0]].vt[0].c[0] += (min_dist_fin[0]/2);
  link_join->cp[cn[0]].vt[0].c[1] += (min_dist_fin[0]/2);
  
  //Now we copy all the vertices in between...

  for (y = min_dist_vert[0]+1; y < link[0]->cp[cn[0]].nv - dflag; y++) {//to the end of the component...
    vect_copy(&(link_join->cp[cn[0]].vt[y-min_dist_vert[0]]),&(link[0]->cp[cn[0]].vt[y]));
  }
  for (y = 0; y < min_dist_vert[0]; y++) {//back to the beginning of the component...
    vect_copy(&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-min_dist_vert[0]-dflag+y]),&(link[0]->cp[cn[0]].vt[y]));
  }

  //And now we're back to creating the other split vertex, this time by
  //copying the original guy and then moving him -rad/sqrt(2) in the
  // x = y  direction.

  vect_copy (&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag]),&(link[0]->cp[cn[0]].vt[min_dist_vert[0]]));
  link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag].c[0] -= (min_dist_fin[0]/2);
  link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag].c[1] -= (min_dist_fin[0]/2);
  
  // So where I once had ------o-----
  //
  // I now have:                  o
  //                      ----     `----
  //                          `o


  //Now I add in my first bridge from this last split vertex
  //to a split vertex on link[1]. In order to do this,
  //I calculate which way I should draw my bridge. Since I can connect
  //a bridging segment to one of two split vertices on the other knot,
  //I want to connect things so that I get:
                  
  //                    /--------
  //                   /
  //                  o
  //              
  //             o    
  //       -----/      

  //instead of 
  //                    /--------
  //		       /
  //                  o
  //                 /|
  //                o |
  //      ------------|
                   
  //in other words, I'm picking which vertex to bridge with which
  //based on what won't risk unwanted crossings. I do this by picking
  //which bridge will allow for link[1]'s new vertices to join up
  //using the shortest segments.

  link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[0] += (min_dist_fin[0]/2);
  link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[1] -= (min_dist_fin[0]/2);

  sel1_dist = threed_len(&(link[1]->cp[cn[1]].vt[min_dist_vert[1]]),&(link[1]->cp[cn[1]].vt[min_dist_vert[1]+1]));

  link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[0] -= min_dist_fin[0];
  link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[1] += min_dist_fin[0];

  sel2_dist = threed_len(&(link[1]->cp[cn[1]].vt[min_dist_vert[1]]),&(link[1]->cp[cn[1]].vt[min_dist_vert[1]+1]));
  
  link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[0] += (min_dist_fin[0]/2);
  link[1]->cp[cn[1]].vt[min_dist_vert[1]].c[1] -= (min_dist_fin[0]/2);

  //draw in my bridging segments...

  for (y = 1; y < new_segs; y++) {

    vect_copy (&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y]),&(link[0]->cp[cn[0]].vt[min_dist_vert[0]]));
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y].c[2] -= y*average;

    //So now I know that if I draw in the bridges by with variation in the xy-plane,
    //I have a better fit.
    if(sel1_dist < sel2_dist) {
      
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y].c[0] -= (min_dist_fin[0]/2);
      //add the x back in gradually.
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y].c[0] += ((y-1)*min_dist_fin[0])/(new_segs-1); 
      
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y].c[1] -= (min_dist_fin[0]/2);
      
      //Here the yz-plane gives the better fit.
    } else {

      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y].c[0] -= (min_dist_fin[0]/2);
      
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y].c[1] -= (min_dist_fin[0]/2);
      //add the y back in gradually.
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv-dflag+y].c[1] += ((y-1)*min_dist_fin[0])/(new_segs-1); 

    }
  }

  //my first bridge has been created. I'll copy in the first split vertex created
  //from link[1]'s min_dist_vert.

  vect_copy (&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv - dflag + new_segs]),&(link[1]->cp[cn[1]].vt[min_dist_vert[1]]));

  if(sel1_dist < sel2_dist) {
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv - dflag + new_segs].c[0] += (min_dist_fin[0]/2);
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv - dflag + new_segs].c[1] -= (min_dist_fin[0]/2);
  } else {
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv - dflag + new_segs].c[0] -= (min_dist_fin[0]/2);
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv - dflag + new_segs].c[1] += (min_dist_fin[0]/2);
  }

  //now I copy in all the unchanged vertices to the end of the array...

  for (y = min_dist_vert[1]+1; y < link[1]->cp[cn[1]].nv - dflag; y++) {
    vect_copy(&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv + new_segs + y - dflag - min_dist_vert[1]]),&(link[1]->cp[cn[1]].vt[y]));
  }
  
  //and I loop back around

  for (y = 0; y < min_dist_vert[1]; y++) {
    vect_copy(&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv -(2*dflag) + new_segs + link[1]->cp[cn[1]].nv-min_dist_vert[1]+y]),&(link[1]->cp[cn[1]].vt[y]));
  }

  //now I'm copying in the other of link[1]'s split vertices.

  vect_copy (&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs]),&(link[1]->cp[cn[1]].vt[min_dist_vert[1]]));
  
  if(sel1_dist < sel2_dist) {
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs].c[0] -= (min_dist_fin[0]/2);
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs].c[1] += (min_dist_fin[0]/2);
  } else {
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs].c[0] += (min_dist_fin[0]/2);
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs].c[1] -= (min_dist_fin[0]/2);
  }

  //and my second bridge to close off the new link.

  for (y = 1; y < new_segs; y++) {
    vect_copy (&(link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y]),&(link[1]->cp[cn[1]].vt[min_dist_vert[1]]));
    link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y].c[2] += y*average;

    if(sel1_dist < sel2_dist) {
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y].c[0] -= (min_dist_fin[0]/2);
      //add the x back in gradually
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y].c[0] += ((y-1)*min_dist_fin[0])/(new_segs-1);
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y].c[1] += (min_dist_fin[0]/2);

    } else {

      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y].c[0] += (min_dist_fin[0]/2);
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y].c[1] -= (min_dist_fin[0]/2);
      //add the y back in gradually
      link_join->cp[cn[0]].vt[link[0]->cp[cn[0]].nv+link[1]->cp[cn[1]].nv - (2*dflag) + new_segs+y].c[1] += ((y-1)*min_dist_fin[0])/(new_segs-1);

    }

  }
  //AND WE'RE DONE.

  plc_write(outfile,link_join);

  free(new_nv);
  free(new_acyclic);
  free(new_colors);

  free(high_low1);
  free(high_low2);

  plc_free(link[0]);
  plc_free(link[1]);
  free(link);

  plc_free(link_join);
			       
  fclose(infile1);
  fclose(infile2);
  fclose(outfile);

  return 1;

}

/* //main() just does some sanity checking and help line printing. */

/* int main(int argv, char *argc[]) { */
/* //printf("%i", 1); */
/*   int i; */
/*   char d = 'e'; */
  
/*   char rcs_revision[128] = "$Revision: 1.10 $"; */
/*   char revtag[128]; */

/*   if (sscanf(&(rcs_revision[1]),"Revision: %s $",revtag) != 1)  */
/*     {strcpy(revtag,rcs_revision);} */
  

/*   if (argv < 6) { */
/*     printf("%i", argv); */
/*     printf("*****************************************************************************************\n"); */
/*     printf("Welcome to knot_add %s by John Foreman (using libraries written by Jason Cantarella et. al.)\nUniversity of Georgia, jforeman@uga.edu\n",revtag);  */
/*     printf("Syntax for knot_add:\n"); */
/*     printf("./knot_add PATH/input1.vect PATH/input2.vect PATH/output.vect component1 component2\n"); */
/*     printf("Example: ./knot_add ./3_1.vect ./5_1.vect ./3_1#5_1.vect 0 0\n"); */
/*     printf("Note that both components are 0 since the knot has only one component...the zeroeth.\n"); */
/*     printf("*****************************************************************************************\n"); */
/*   }  */
/*   else { */
/*     for(i=0; i<6; i++){ */
/*       printf("%c", argc[i][0]);	 */
/* 	} */
/*   main_func(argv, argc); */
/*   } */

/* } */


