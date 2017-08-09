/*******************************************************************

 orient.c : Reorients an OFF file to provide a "nice" view down the x-axis.
            If multiple files are given, splines the orientations to provide a 
            smooth reorientation of the curve. 

	   ************************************************************/

#include<gsl/gsl_eigen.h>
#include<gsl/gsl_matrix.h>
#include<plsurf.h>
#include<argtable2.h>
#include<math.h>

/**********************************************************************/

struct arg_lit  *verbose;
struct arg_file *off_file;
struct arg_lit  *help;
struct arg_int  *axis;
struct arg_int  *spinframes; /* Adds a slow spin (of this many frames)
				to the last file */

struct arg_end *end;
struct arg_end *helpend;

surface       surf;
FILE         *infile_fptr,*outfile_fptr;
int           VERBOSE = 0;

typedef struct quattype {

  double x,y,z,w;

} quat;

void recenter_surf(surface *surf)
/* Translate to put center of mass at the origin. */
{
  int i;

  if (VERBOSE > 5) {
    
    printf("Computing center of mass for %d vertex surface...\n",surf->verts);
    
  }
  
  plc_vector com;
  com = (plc_vector)(vertex_center_of_mass(surf));

  for(i=0;i<surf->verts;i++) {

    surf->vert_buf[i] = plc_vect_diff(surf->vert_buf[i],com);

  }

}

void reframe_surf(surface *surf,plc_vector rows[3],int axis) 
{ 
  int i,j;
  plc_vector tempV;

  for(i=0;i<surf->verts;i++) {
    
    tempV = surf->vert_buf[i];
    
    for(j=0;j<3;j++) {

      surf->vert_buf[i].c[j] = plc_dot_prod(rows[(axis+j)%3],tempV);

    }

  }
}


double ReciprocalSqrt( float x ) { 
  /*long i; 
    float y, r; 
 
    y = x * 0.5f; 
    i = *(long *)( &x ); 
    i = 0x5f3759df - ( i >> 1 ); 
    r = *(float *)( &i ); 
    r = r * ( 1.5f - r * r * y ); */ 
  return 1/sqrt(x); 
} 

quat matrix_to_quaternion3( plc_vector rows[3]) {

  quat q;

  q.w = sqrt(rows[0].c[0] + rows[1].c[1] + rows[2].c[2] + 1.0)/2.0;
  q.x = (rows[2].c[1] - rows[1].c[2]) / (4.0 * q.w) ;
  q.y = (rows[0].c[2] - rows[2].c[0]) / (4.0 * q.w) ;
  q.z = (rows[1].c[0] - rows[0].c[1]) / (4.0 * q.w) ;

  return q;
}

quat matrix_to_quaternion2( plc_vector rows[3] ) { 
 
  quat q;

  if ( rows[0].c[0] /* m[0 * 4 + 0] */ + rows[1].c[1] /* m[1 * 4 + 1] */ + rows[2].c[2] /*m[2 * 4 + 2]*/ > 0.0f ) { 
    
    float t = + rows[0].c[0] /*m[0 * 4 + 0]*/ + rows[1].c[1] /*m[1 * 4 + 1]*/+ rows[2].c[2] /*m[2 * 4 + 2]*/ + 1.0f; 
    float s = ReciprocalSqrt( t ) * 0.5f; 
    
    q.w = s * t; 
    q.z = ( rows[0].c[1] /*m[0 * 4 + 1]*/ - rows[1].c[0] /*m[1 * 4 + 0]*/ ) * s; 
    q.y = ( rows[2].c[0] /*m[2 * 4 + 0]*/ - rows[0].c[2] /*m[0 * 4 + 2]*/ ) * s; 
    q.x = ( rows[1].c[2] /*m[1 * 4 + 2]*/ - rows[2].c[1] /*m[2 * 4 + 1]*/ ) * s; 
    
  } else if ( rows[0].c[0] /*m[0 * 4 + 0]*/ > rows[1].c[1] /*m[1 * 4 + 1]*/ && rows[0].c[0] /*m[0 * 4 + 0]*/ > rows[2].c[2] /*m[2 * 4 + 2]*/ ) { 
    
    float t = + rows[0].c[0] /*m[0 * 4 + 0]*/ - rows[1].c[1] /*m[1 * 4 + 1]*/ - rows[2].c[2] /*m[2 * 4 + 2]*/ + 1.0f; 
    float s = ReciprocalSqrt( t ) * 0.5f; 
    
    q.x = s * t; 
    q.y = ( rows[0].c[1] /*m[0 * 4 + 1]*/ + rows[1].c[0] /*m[1 * 4 + 0]*/ ) * s; 
    q.z = ( rows[2].c[0] /*m[2 * 4 + 0]*/ + rows[0].c[2] /*m[0 * 4 + 2]*/ ) * s; 
    q.w = ( rows[1].c[2] /*m[1 * 4 + 2]*/ - rows[2].c[1] /*m[2 * 4 + 1]*/ ) * s; 
    
  } else if ( rows[1].c[1] /*m[1 * 4 + 1]*/ > rows[2].c[2] /*m[2 * 4 + 2]*/ ) { 
    
    float t = -rows[0].c[0] /*m[0 * 4 + 0]*/ + rows[1].c[1] /*m[1 * 4 + 1]*/ - rows[2].c[2] /*m[2 * 4 + 2]*/ + 1.0f; 
    float s = ReciprocalSqrt( t ) * 0.5f; 
    
    q.y = s * t; 
    q.x = ( rows[0].c[1] /*m[0 * 4 + 1]*/ + rows[1].c[0] /*m[1 * 4 + 0]*/ ) * s; 
    q.w = ( rows[2].c[0] /*m[2 * 4 + 0]*/ - rows[0].c[2] /*m[0 * 4 + 2]*/ ) * s; 
    q.z = ( rows[1].c[2] /*m[1 * 4 + 2]*/ + rows[2].c[1] /*m[2 * 4 + 1]*/ ) * s; 
    
  } else { 
    
    float t = - rows[0].c[0] /*m[0 * 4 + 0]*/ - rows[1].c[1]/*m[1 * 4 + 1]*/ + rows[2].c[2]/*m[2 * 4 + 2]*/ + 1.0f; 
    float s = ReciprocalSqrt( t ) * 0.5f; 
    
    q.z = s * t; 
    q.w = ( rows[0].c[1] /*m[0 * 4 + 1]*/ - rows[1].c[0] /*m[1 * 4 + 0]*/ ) * s; 
    q.x = ( rows[2].c[0] /*m[2 * 4 + 0]*/ + rows[0].c[2] /*m[0 * 4 + 2]*/ ) * s; 
    q.y = ( rows[1].c[2] /*m[1 * 4 + 2]*/ + rows[2].c[1] /*m[2 * 4 + 1]*/ ) * s; 
    
  } 
 
  return q;

} 

quat matrix_to_quaternion(plc_vector rows[3])
/* Convert matrix to quaternion representation */
{
  quat q;
  double m00, m01, m02, m10, m11, m12, m20, m21, m22;

  m00 = rows[0].c[0];  m10 = rows[1].c[0];  m20 = rows[2].c[0];
  m01 = rows[0].c[1];  m11 = rows[1].c[1];  m21 = rows[2].c[1];
  m02 = rows[0].c[2];  m12 = rows[1].c[2];  m22 = rows[2].c[2];

  double tr = m00 + m11 + m22;
  
  if (tr > 0) { 
    double S = sqrt(tr+1.0) * 2; // S=4*qw 
    q.w = 0.25 * S;
    q.x = (m21 - m12) / S;
    q.y = (m02 - m20) / S; 
    q.z = (m10 - m01) / S; 
  } else if ((m00 > m11) &&  (m00 > m22)) { 
    double S = sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx 
    q.w = (m21 - m12) / S;
    q.x = 0.25 * S;
    q.y = (m01 + m10) / S; 
    q.z = (m02 + m20) / S; 
  } else if (m11 > m22) { 
    double S = sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
    q.w = (m02 - m20) / S;
    q.x = (m01 + m10) / S; 
    q.y = 0.25 * S;
    q.z = (m12 + m21) / S; 
  } else { 
    double S = sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
    q.w = (m10 - m01) / S;
    q.x = (m02 + m20) / S;
    q.y = (m12 + m21) / S;
    q.z = 0.25 * S;
  }

  return q;
}

quat slerp(quat qa, quat qb, double t) {
  // quaternion to return
  quat qm;
  // Calculate angle between them.
  double cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;
  // if qa=qb or qa=-qb then theta = 0 and we can return qa
  if (abs(cosHalfTheta) >= 1.0){
    qm.w = qa.w;qm.x = qa.x;qm.y = qa.y;qm.z = qa.z;
    return qm;
  }
  // Calculate temporary values.
  double halfTheta = acos(cosHalfTheta);
  double sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);
  // if theta = 180 degrees then result is not fully defined
  // we could rotate around any axis normal to qa or qb
  if (fabs(sinHalfTheta) < 0.001){ // fabs is floating point absolute
    qm.w = (qa.w * 0.5 + qb.w * 0.5);
    qm.x = (qa.x * 0.5 + qb.x * 0.5);
    qm.y = (qa.y * 0.5 + qb.y * 0.5);
    qm.z = (qa.z * 0.5 + qb.z * 0.5);
    return qm;
  }
  double ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
  double ratioB = sin(t * halfTheta) / sinHalfTheta; 
  //calculate Quaternion.
  qm.w = (qa.w * ratioA + qb.w * ratioB);
  qm.x = (qa.x * ratioA + qb.x * ratioB);
  qm.y = (qa.y * ratioA + qb.y * ratioB);
  qm.z = (qa.z * ratioA + qb.z * ratioB);
  return qm;
}

plc_vector *quaternion_to_matrix(quat q)
/* Convert quaternion to matrix */
{
  plc_vector *rows;
  rows = calloc(3,sizeof(plc_vector));

  double x2 = q.x + q.x; 
  double y2 = q.y + q.y; 
  double z2 = q.z + q.z; 
  { 
    double xx2 = q.x * x2; 
    double yy2 = q.y * y2; 
    double zz2 = q.z * z2; 
    
    rows[0].c[0] = 1.0 - yy2 - zz2; 
    rows[1].c[1] = 1.0 - xx2 - zz2; 
    rows[2].c[2] = 1.0 - xx2 - yy2; 
  } 
  { 
    double yz2 = q.y * z2; 
    double wx2 = q.w * x2; 
    
    rows[2].c[1] = yz2 - wx2; 
    rows[1].c[2] = yz2 + wx2; 
  } 
  { 
    double xy2 = q.x * y2; 
    double wz2 = q.w * z2; 
    
    rows[1].c[0] = xy2 - wz2; 
    rows[0].c[1] = xy2 + wz2; 
  } 
  { 
    double xz2 = q.x * z2; 
    double wy2 = q.w * y2; 
    
    rows[0].c[2] = xz2 - wy2; 
    rows[2].c[0] = xz2 + wy2; 
  } 

return rows;
}

plc_vector *quaternion_to_matrix2(quat q)
/* Convert quaternion to matrix */
{
  plc_vector *rows;
  rows = calloc(3,sizeof(plc_vector));

  double xx,yy,zz,ww,xy,xz,xw,yz,yw,zw;

  xx = q.x * q.x;
  yy = q.y * q.y;
  zz = q.z * q.z;
  ww = q.w * q.w;

  xy = q.x * q.y;
  xz = q.x * q.z;
  xw = q.x * q.w;

  yz = q.y * q.z;
  yw = q.y * q.w;

  zw = q.z * q.w;

  rows[0].c[0]  = 1 - 2 * ( yy + zz );
  rows[0].c[1]  =     2 * ( xy - zw );
  rows[0].c[2] =     2 * ( xz + yw );
  
  rows[1].c[0]  =     2 * ( xy + zw );
  rows[1].c[1]  = 1 - 2 * ( xx + zz );
  rows[1].c[2]  =     2 * ( yz - xw );
  
  rows[2].c[0]  =     2 * ( xz - yw );
  rows[2].c[1]  =     2 * ( yz + xw );
  rows[2].c[2] = 1 - 2 * ( xx + yy );

  return rows;
}

// typedef struct {double x,y,z,w;} quat; 
typedef double HMatrix[4][4]; 
#define X 0 
#define Y 1 
#define Z 2 
#define W 3 
/* Return quaternion product qL * qR. */ 
quat Qt_Mul(quat qL, quat qR) 
{ 
quat qq; 
qq.w = qL.w*qR.w - qL.x*qR.x - qL.y*qR.y - qL.z*qR.z; 
qq.x = qL.w*qR.x + qL.x*qR.w + qL.y*qR.z - qL.z*qR.y; 
qq.y = qL.w*qR.y + qL.y*qR.w + qL.z*qR.x - qL.x*qR.z; 
qq.z = qL.w*qR.z + qL.z*qR.w + qL.x*qR.y - qL.y*qR.x; 
return (qq); 
} 
/* Return norm of quaternion, the sum of the squares of the components. */ 
#define Qt_Norm(q) ((q).x*(q).x + (q).y*(q).y + (q).z*(q).z + (q).w*(q).w) 
/* Construct rotation matrix from (possibly non-unit) quaternion. 
 * Assumes matrix is used to multiply column vector on the left: 
 * vnew = mat vold.  Works correctly for right-handed coordinate system 
 * and right-handed rotations. */ 
void Qt_ToMatrix(quat q, HMatrix mat) 
{ 
double Nq = Qt_Norm(q); 
double s = (Nq > 0.0) ? (2.0 / Nq) : 0.0; 
double xs = q.x*s, ys = q.y*s, zs = q.z*s; 
double wx = q.w*xs, wy = q.w*ys, wz = q.w*zs; 
double xx = q.x*xs, xy = q.x*ys, xz = q.x*zs; 
double yy = q.y*ys, yz = q.y*zs, zz = q.z*zs; 
mat[X][X] = 1.0 - (yy + zz); mat[Y][X] = xy + wz; mat[Z][X] = xz - wy; 
mat[X][Y] = xy - wz; mat[Y][Y] = 1.0 - (xx + zz); mat[Z][Y] = yz + wx; 
mat[X][Z] = xz + wy; mat[Y][Z] = yz - wx; mat[Z][Z] = 1.0 - (xx + yy); 
mat[X][W] = mat[Y][W] = mat[Z][W] = 0.0; 
mat[W][X] = mat[W][Y] = mat[W][Z] = 0.0; 
mat[W][W] = 1.0; 
}

/* Construct a unit quaternion from rotation matrix.  Assumes matrix is 
 * used to multiply column vector on the left: vnew = mat vold.  Works 
 * correctly for right-handed coordinate system and right-handed rotations. 
 * Translation and perspective components ignored. */ 

quat Qt_FromMatrix(HMatrix mat) 
{ 
/* This algorithm avoids near-zero divides by looking for a large component 
 * first w, then x, y, or z.  When the trace is greater than zero, 
 * |w| is greater than 1/2, which is as small as a largest component can be. 
 * Otherwise, the largest diagonal entry corresponds to the largest of |x|, 
 * |y|, or |z|, one of which must be larger than |w|, and at least 1/2. */ 
quat qu = {0,0,0,0}; 
double tr, s; 
tr = mat[X][X] + mat[Y][Y]+ mat[Z][Z]; 
if (tr >= 0.0) { 
s = sqrt(tr + mat[W][W]); 
qu.w = s*0.5; 
s = 0.5 / s; 
qu.x = (mat[Z][Y] - mat[Y][Z]) * s; 
qu.y = (mat[X][Z] - mat[Z][X]) * s; 
qu.z = (mat[Y][X] - mat[X][Y]) * s; 
} else { 
int h = X; 
if (mat[Y][Y] > mat[X][X]) h = Y; 
if (mat[Z][Z] > mat[h][h]) h = Z; 
switch (h) { 
#define caseMacro(i,j,k,I,J,K) \
case I:\
s = sqrt( (mat[I][I] - (mat[J][J]+mat[K][K])) + mat[W][W] );\
qu.i = s*0.5;\
s = 0.5 / s;\
qu.j = (mat[I][J] + mat[J][I]) * s;\
qu.k = (mat[K][I] + mat[I][K]) * s;\
qu.w = (mat[K][J] - mat[J][K]) * s;\
break
caseMacro(x,y,z,X,Y,Z); 
caseMacro(y,z,x,Y,Z,X); 
caseMacro(z,x,y,Z,X,Y); 

} 
} 
if (mat[W][W] != 1.0) { 
s = 1.0/sqrt(mat[W][W]); 
qu.w *= s; qu.x *= s; qu.y *= s; qu.z *= s; 
} 
return (qu); 
} 

quat matrix_to_quaternion4(plc_vector frameA[3])
{
  HMatrix H = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}};
  int i,j;

  for(i=0;i<3;i++) { 
    for(j=0;j<3;j++) {
      H[i][j] = frameA[i].c[j];
    }
  }

  return Qt_FromMatrix(H);
}

plc_vector *quaternion_to_matrix4(quat q)
{
  HMatrix H;
  plc_vector *frameA = calloc(3,sizeof(plc_vector));
  Qt_ToMatrix(q, H);

  int i,j;
  
  for(i=0;i<3;i++) { 
    for(j=0;j<3;j++) {
      frameA[i].c[j] = H[i][j];
    }
  }
  
  return frameA;
}
  

plc_vector *minterp(plc_vector frameA[3],plc_vector frameB[3],double t)
{
  plc_vector *rows;
  quat qA,qB,qS;

  if (fabs(t) < 1e-10) { 

    rows = calloc(3,sizeof(plc_vector));
    int i;
    for(i=0;i<3;i++) { rows[i] = frameA[i]; }
    
  } else {

    qA = matrix_to_quaternion4(frameA);
    qB = matrix_to_quaternion4(frameB);
    qS = slerp(qA,qB,t);
    rows = quaternion_to_matrix4(qS);
    
  }

  return rows;
} 

plc_vector *inertialframe_surf(surface *surf)
/* Compute inertial axes of surface */
{
  gsl_matrix *I = gsl_matrix_calloc(3,3); 
  int i,j,k;
  double temp;
  plc_vector *rows;

  rows = calloc(3,sizeof(plc_vector));
 
  if (VERBOSE > 5) {

    printf("Computing inertial tensor...");

  }
  
  /* Now we compute the inertial tensor */


  for (i=0;i<surf->verts;i++) {

    for (j=0;j<3;j++) {

      for(k=0;k<3;k++) {

	temp = gsl_matrix_get(I,j,k);
	temp += surf->vert_buf[i].c[j] * surf->vert_buf[i].c[k];
	gsl_matrix_set(I,j,k,temp);

      }

    }

  }
	
  gsl_eigen_symmv_workspace *wks = gsl_eigen_symmv_alloc(3);
  gsl_matrix *Evec = gsl_matrix_calloc(3,3);
  gsl_vector *Eval = gsl_vector_calloc(3);

  if (VERBOSE > 5) {

    printf("done.\n");
    printf("Inertia matrix %4.2g %4.2g %4.2g\n"
	   "               %4.2g %4.2g %4.2g\n"
	   "               %4.2g %4.2g %4.2g\n\n",
	   gsl_matrix_get(I,0,0),gsl_matrix_get(I,0,1),gsl_matrix_get(I,0,2),
	   gsl_matrix_get(I,1,0),gsl_matrix_get(I,1,1),gsl_matrix_get(I,1,2),
	   gsl_matrix_get(I,2,0),gsl_matrix_get(I,2,1),gsl_matrix_get(I,2,2));

    printf("Diagonalizing I matrix...");

  }

  gsl_eigen_symmv(I,Eval,Evec,wks);
  gsl_eigen_symmv_sort(Eval,Evec,GSL_EIGEN_SORT_ABS_DESC);
  
  if (VERBOSE > 5) {

    printf("done.\n");
    printf("Axes of rotation are\n %g %g %g \n %g %g %g \n %g %g %g \n\n",
	   gsl_matrix_get(Evec,0,0),gsl_matrix_get(Evec,0,1),gsl_matrix_get(Evec,0,2),
	   gsl_matrix_get(Evec,1,0),gsl_matrix_get(Evec,1,1),gsl_matrix_get(Evec,1,2),
	   gsl_matrix_get(Evec,2,0),gsl_matrix_get(Evec,2,1),gsl_matrix_get(Evec,2,2));
    printf("Eigenvalues are %g, %g, %g.\n",
	   gsl_vector_get(Eval,0),gsl_vector_get(Eval,1),gsl_vector_get(Eval,2));

    printf("Rotating to make these the x, y, and z axes...\n");

  }

  /* We have now computed the principal axes of rotation and sorted
     them.  We want to change coordinates so that the eigenvectors
     become the x,y, and z axes, which means multiplying by the
     inverse of Evec. Luckily, Evec is an orthgonal matrix, so we can
     just take the transpose of Evec. */
  
  gsl_matrix_transpose(Evec);

  for(i=0;i<3;i++) {

    for(j=0;j<3;j++) {
      
      rows[i].c[j] = gsl_matrix_get(Evec,j,i);
      
    }

  }

  /* We may now have a problem. The eigenvectors are chosen at random, and so 
     we may actually have a mirror/rotation matrix rather than a rotation matrix. The 
     best solution is to check the determinant and flip if needed. */

  if (plc_dot_prod(rows[0],plc_cross_prod(rows[1],rows[2])) < 0) {

    rows[1] = plc_scale_vect(-1.0,rows[1]);

  }

  gsl_matrix_free(Evec);
  gsl_matrix_free(I);
  gsl_vector_free(Eval);

  return rows;
 
}

surface loadsurf(const char *filename) 
{

  surface surf;
  FILE *infile_fptr;
 
  /* We begin by loading the current surface. */
    
  infile_fptr = fopen(filename,"r");
  
  if (infile_fptr == NULL) {
    
      fprintf(stderr,"orient: Couldn't open file %s.\n",filename);
      exit(1);
      
  }
  
  if (!load_surf_from_OFF(&surf,infile_fptr)) {
    
    fprintf(stderr,"orient: Couldn't read OFF file %s.\n",filename);
    exit(1);  /* Try the next file */
    
  }
  
  fclose(infile_fptr);
  return surf;
}

void transpose(plc_vector rows[3]) {

  plc_vector cr[3];
  int i;

  for(i=0;i<3;i++) {cr[i] = rows[i];}

  rows[0].c[1] = cr[1].c[0];
  rows[0].c[2] = cr[2].c[0];
  rows[1].c[0] = cr[0].c[1];
  rows[1].c[2] = cr[2].c[1];
  rows[2].c[0] = cr[0].c[2];
  rows[2].c[1] = cr[1].c[2];

}
  

/****************************** Main procedure ********************************/
  
int main(int argc,char *argv[])
{
  int            infilenum,nerrors;
  int            daxis = 0;
 
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     off_file = arg_filen(NULL,NULL,"<file>",1,10000,"input files"),
     axis = arg_int0("a","axis","<0-2>","which axis of inertia to display"),
     spinframes = arg_int0("s","spinframes","<n>","adds spin to last model over n frames"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("orient: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"orient");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("orient reorients OFF files to present a best view down the x axis.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (verbose -> count > 0) {

    VERBOSE = 10;

  }

  if (axis -> count > 0) {

    daxis = axis->ival[0];

  }

  if (daxis < 0 || daxis > 2) {

    printf("reorient: Axis must be 0, 1, or 2.\n");
    exit(0);

  }

  plc_vector *frameStart, *frameEnd;
  surface surfStart, surfEnd;
  /* We need to find the frame for the starting and ending files. */

  surfStart = loadsurf(off_file->filename[0]);
  surfEnd = loadsurf(off_file->filename[off_file->count-1]);

  recenter_surf(&surfStart);
  recenter_surf(&surfEnd);

  frameStart = inertialframe_surf(&surfStart);
  frameEnd = inertialframe_surf(&surfEnd);
  
  /* Now we have parsed the arguments and are ready to work. */
 
  double t;
  plc_vector *frame;

  for(t=0,infilenum = 0;infilenum < off_file->count;infilenum++,t = ((double)(infilenum)/(double)(off_file->count))) {

    surf = loadsurf(off_file->filename[infilenum]);
    recenter_surf(&surf);
    
    frame = minterp(frameStart,frameEnd,t);

    if (VERBOSE > 5) {

      printf("\nChecking interpolated frame...\n");

      printf(" %g %g %g \n %g %g %g \n %g %g %g\n\n",
	     plc_M_clist(frame[0]),
	     plc_M_clist(frame[1]),
	     plc_M_clist(frame[2]));
    }

    reframe_surf(&surf,frame,daxis);
    free(frame);

    plc_vector *frame_check;
    frame_check = inertialframe_surf(&surf);

    /* At this point, the inertial frame ought to be aligned with the axes */

    if (VERBOSE > 5) {

      printf("\nChecking inertial frame after rotation is applied...\n");

      printf(" %g %g %g \n %g %g %g \n %g %g %g\n\n",
	     plc_M_clist(frame_check[0]),
	     plc_M_clist(frame_check[1]),
	     plc_M_clist(frame_check[2]));
    }

    /* Now overwrite the original file. */

    outfile_fptr = fopen(off_file->filename[infilenum],"w");
    write_surf_to_OFF(&surf,outfile_fptr);
    fclose(outfile_fptr);

    /* Now do housekeeping. */

    kill_surface(&surf);

  }

  infilenum--;

  if (spinframes->count > 0) {

    printf("Writing %d spin frames...",spinframes->ival[0]);

    int i;
    
    transpose(frameEnd);

    for(i=0,t=0;i<spinframes->ival[0];i++) {

      double x,PI=3.1415926;
      x = (double)(i)/(double)(spinframes->ival[0]-1);
      t = 2*PI*((cos(PI + PI*x) + 1)/2.0); /* from 0 to 1, with acceleration */

      plc_vector newframe[3];

      newframe[0] = plc_vlincomb(cos(t),frameEnd[0],sin(t),frameEnd[2]);
      newframe[1] = frameEnd[1];
      newframe[2] = plc_vlincomb(cos(t+PI/2.0),frameEnd[0],
				 sin(t+PI/2.0),frameEnd[2]);

      transpose(newframe);

      surf = loadsurf(off_file->filename[infilenum]);
      recenter_surf(&surf);
      reframe_surf(&surf,newframe,2);
      
      /* Now write the spin file. */

      char filename[1024];
      FILE *outfile;

      sprintf(filename,"spin.%07d.off",i);
      outfile = fopen(filename,"w");
      write_surf_to_OFF(&surf,outfile);
      fclose(outfile);
      
      /* Now do housekeeping. */
      
      kill_surface(&surf);

    }

    printf("done.\n");
 
  }

  exit(0);

}




