/*
 *  anymodel.c
 *  
 *
 *  Created by cantarel on Thu Apr 17 2003.
 */

#include"writhecomp.h"

NUMPTR  np_plusequal(NUMPTR a,NUMPTR b)

     /* Procedure computes a += b, with memory management. */
     /* It is explicitly safe to call 

              a = np_plusequal(a,b);

     */

{
  NUMPTR temp;

  temp = (*np_add)(a,b);
  (*kill_np)(a);
  
  return temp;
}

NUMPTR  np_minusequal(NUMPTR a,NUMPTR b)
    
     /* Procedure computes a -= b, with memory management. */
     /* It is explicitly safe to call 

              a = np_minusequal(a,b);

     */

{ 
  NUMPTR temp,temp2;

  temp2 = (*np_mul)(np_MINUS,b);
  temp  = (*np_add)(a,temp2);

  (*kill_np)(a);
  (*kill_np)(temp2);
  
  return temp;
}

NUMPTR *vector_to_nv(plc_vector a)

     /* Procedure converts the vector a to a NUMVECTOR. */

{
  NUMPTR *A;
 
  A = calloc(3,sizeof(NUMPTR));

  A[0] = double_to_np(a.c[0]);
  A[1] = double_to_np(a.c[1]);
  A[2] = double_to_np(a.c[2]);

  return A;
}

NUMPTR *nv_diff(NUMVECTOR a,NUMVECTOR b)

/* Procedure subtracts two numvectors. */

{
    NUMPTR   *ret;
    NUMPTR    tmp;
    
    int i;

    ret = calloc(3,sizeof(NUMPTR));
  
    for(i=0;i<3;i++) {

        tmp    = (*np_mul)(np_MINUS,b[i]);
        ret[i] = (*np_add)(a[i],tmp);
        (*kill_np)(tmp);

    }

    return ret;
}

NUMPTR *nv_cross(NUMVECTOR a, NUMVECTOR b)

/* Procedure computes the cross product of a and b. */

{
    NUMPTR    *ret;
    NUMPTR    tmp[3];

    ret = calloc(3,sizeof(NUMPTR));

    /* Compute first component */
    
    tmp[0] = (*np_mul)(a[1],b[2]);
    tmp[1] = (*np_mul)(a[2],b[1]);
    tmp[2] = (*np_mul)(np_MINUS,tmp[1]);
    
    ret[0] = (*np_add)(tmp[0],tmp[2]);

    (*kill_np)(tmp[0]); (*kill_np)(tmp[1]); (*kill_np)(tmp[2]);

    /* Compute second component */
    
    tmp[0] = (*np_mul)(a[2],b[0]);
    tmp[1] = (*np_mul)(a[0],b[2]);
    tmp[2] = (*np_mul)(np_MINUS,tmp[1]);

    ret[1] = (*np_add)(tmp[0],tmp[2]);

    (*kill_np)(tmp[0]); (*kill_np)(tmp[1]); (*kill_np)(tmp[2]);

    /* Compute third component */

    tmp[0] = (*np_mul)(a[0],b[1]);
    tmp[1] = (*np_mul)(a[1],b[0]);
    tmp[2] = (*np_mul)(np_MINUS,tmp[1]);

    ret[2] = (*np_add)(tmp[0],tmp[2]);

    (*kill_np)(tmp[0]); (*kill_np)(tmp[1]); (*kill_np)(tmp[2]);

    return ret;
}

NUMPTR nv_dot(NUMVECTOR a, NUMVECTOR b)

/* Procedure computes the dot product of two vectors. */

{
    NUMPTR ret;
    NUMPTR tmp[4];
    int    i;

    tmp[0] = (*np_mul)(a[0],b[0]);
    tmp[1] = (*np_mul)(a[1],b[1]);
    tmp[2] = (*np_mul)(a[2],b[2]);

    tmp[3] = (*np_add)(tmp[0],tmp[1]);

    ret = (*np_add)(tmp[3],tmp[2]);

    for(i=0;i<4;i++) (*kill_np)(tmp[i]);

    return ret;
}

void nv_normalize(NUMVECTOR a)

/* Procedure computes the normalization of a. */

{
    NUMVECTOR tmp;
    NUMPTR    norm,dt;
    int       i;

    dt   = nv_dot(a,a);
    norm = (*np_sqrt)(dt);
    (*kill_np)(dt);

    /* We need to do a little magic here in order to change the
        values of a "in place".

        Notice that the values in a right now are POINTERS, which
        means we can free and discard them as we go along, as long
        as we wait until we're done with them. */
    
    for(i=0;i<3;i++) {

        tmp[i] = (*np_div)(a[i],norm);
        (*kill_np)(a[i]);
        a[i] = tmp[i];

    }

    (*kill_np)(norm);

    /* Notice that, in contrast to the other procedures around here,
        we DONT free tmp at the end of the procedure. This isn't a bug,
        since we're returning those pointers in a itself. */
    
}

void kill_nv(NUMVECTOR a)

     /* Procedure frees the memory associated with a by calling kill_np. */

{
  int i;

  for(i=0;i<3;i++) {

    (*kill_np)(a[i]);

  }

}

NUMPTR np_dihedral_sum(NUMVECTOR X_one,NUMVECTOR X_two,
                       NUMVECTOR Y_one,NUMVECTOR Y_two)
/*

 Procedure computes 4 of the dihedral angles of the tetrahedron
 specified by the give vectors and returns their sum. These are all
 dihedrals except those along the edges X_1X_2 and Y_1Y_2.

 We have first must check and make sure that the vectors are not
 coplanar. This is the only degenerate case.

 We return the dihedral sum - twopi, in an effort to support the
 writhe calculation which is our calling procedure. We also assign
 a sign to the crossing based on the determinant of the matrix 

 [ X_two - X_one, Y_two - Y_one, Y_one - X_one ].

 */
{
    NUMPTR *Y_one_X_one, *Y_two_X_one, *Y_two_Y_one,
    *X_two_X_one, *X_two_Y_one;

    NUMPTR *X_one_Y_one_X_two, *X_one_Y_one_Y_two,
        *X_one_Y_two_X_two, *Y_one_Y_two_X_two;

    NUMPTR cosines[4],angles[4],angle_sum,denom,denomsqr;
    NUMPTR tmp[2];
    NUMPTR absdenom, sign, two, twopi, signed_angle_sum;

    int n;

    /* First we compute the located vectors along the edges of the
        tetrahedron that we will need to use. */

    Y_one_X_one = nv_diff(Y_one,X_one);
    Y_two_X_one = nv_diff(Y_two,X_one);
    Y_two_Y_one = nv_diff(Y_two,Y_one);
    X_two_X_one = nv_diff(X_two,X_one);
    X_two_Y_one = nv_diff(X_two,Y_one);

    /* We must now worry that the vectors may be coplanar. We detect this
        by computing the normal to the plane determined by X_One, X_Two,
        and Y_One, and computing the dot product of this normal vector with
        the vector Y_Two - Y_One. If these vectors are perpendicular, the
        tetrahedron is degenerate, and we return appropriate values.  */

    X_one_Y_one_X_two = nv_cross(Y_one_X_one,X_two_X_one);

    /* There is another slightly tricky issue here: How do we know which */
    /* denominators are "too small" in our arithmetic model? We use the  */
    /* constant NUMPTR of np_EPS, which is set by the initializer.      */

    denom    = nv_dot(X_one_Y_one_X_two,Y_two_Y_one);
    denomsqr = np_mul(denom,denom);
    absdenom = np_sqrt(denomsqr);
    
    if(np_leq(absdenom,np_EPS)) {

        /* The vectors are coplanar. We free all the scratch variables */
        /* and return a very good zero. */

      kill_nv(Y_one_X_one); free(Y_one_X_one);
      kill_nv(Y_two_X_one); free(Y_two_X_one);
      kill_nv(Y_two_Y_one); free(Y_two_Y_one);
      kill_nv(X_two_X_one); free(X_two_X_one);
      kill_nv(X_two_Y_one); free(X_two_Y_one);
 
      kill_nv(X_one_Y_one_X_two); free(X_one_Y_one_X_two);

      kill_np(denomsqr);
      kill_np(denom);
      kill_np(absdenom);
      
      tmp[0] = double_to_np(0.0);
      return(tmp[0]);
      
    } else if (np_leq(denom,np_ZERO)) {
      
      /* We now set the sign of the crossing accordingly. */
      
      sign = double_to_np(1.0);
      
    } else {

      sign = double_to_np(-1.0);

    }

    /* Next, we take their cross products to generate normal vector for
        the sides of the tetrahedron.  Before doing so, though, we should
        beware of the many problems which arise in dealing with these cross
        products when they are too small. */

    X_one_Y_one_Y_two = nv_cross(Y_two_X_one,Y_one_X_one);
    X_one_Y_two_X_two = nv_cross(X_two_X_one,Y_two_X_one);
    Y_one_Y_two_X_two = nv_cross(Y_two_Y_one,X_two_Y_one);

    /* Having taken these cross products, the pairwise differences are
        used up, and we can discard them. */

    kill_nv(Y_one_X_one); free(Y_one_X_one);
    kill_nv(Y_two_X_one); free(Y_two_X_one);
    kill_nv(Y_two_Y_one); free(Y_two_Y_one);
    kill_nv(X_two_X_one); free(X_two_X_one);
    kill_nv(X_two_Y_one); free(X_two_Y_one);
    
    /* Now we vectornormalize these vectors in preparation to compute angles. */

    nv_normalize(X_one_Y_one_X_two);
    nv_normalize(X_one_Y_one_Y_two);
    nv_normalize(X_one_Y_two_X_two);
    nv_normalize(Y_one_Y_two_X_two);

    /* Last, we take dot products to compute cosines... */

    cosines[0] = nv_dot(X_one_Y_one_X_two,X_one_Y_one_Y_two);
    cosines[1] = nv_dot(X_one_Y_one_Y_two,X_one_Y_two_X_two);
    cosines[2] = nv_dot(Y_one_Y_two_X_two,X_one_Y_two_X_two);
    cosines[3] = nv_dot(X_one_Y_one_X_two,Y_one_Y_two_X_two);

    /* We have now finished with the cross products, so let's
        go ahead and free them, too. */

    kill_nv(X_one_Y_one_X_two); 
    kill_nv(X_one_Y_one_Y_two);
    kill_nv(X_one_Y_two_X_two);
    kill_nv(Y_one_Y_two_X_two);

    free(X_one_Y_one_X_two); 
    free(X_one_Y_one_Y_two);
    free(X_one_Y_two_X_two);
    free(Y_one_Y_two_X_two);
    
    /* We are now prepared to compute the angles involved. */
    /* A subtle issue arises: we may be slightly over "1.0" */
    /* when we take these cosines. The best thing to do is  */
    /* check the results against np_MINUS, and if they fail, */
    /* to set them accordingly. */


    for(n=0;n < 4;n++) {

        angles[n] = np_acos(cosines[n]);
        kill_np(cosines[n]);

    }

    tmp[0] = np_add(angles[0],angles[1]);
    tmp[1] = np_add(angles[2],angles[3]);     
    angle_sum = np_add(tmp[0],tmp[1]);

    kill_np(tmp[0]); kill_np(tmp[1]); 

    for(n=0;n < 4;n++) {

      kill_np(angles[n]);

    }

    /* We now subtract 2pi, and multiply by "sign" to set the correct */
    /* sign for the crossing. */

    two = double_to_np(2.0);
    twopi = np_mul(two,np_PI);
    angle_sum = np_minusequal(angle_sum,twopi);
    signed_angle_sum = np_mul(angle_sum,sign);

    kill_np(two); kill_np(twopi); kill_np(angle_sum); kill_np(denom); kill_np(sign);
    kill_np(absdenom); kill_np(denomsqr); 

    return( signed_angle_sum );
    
}

NUMPTR np_writhe(plCurve *writhe_pline,int ii,int jj)

     /* Procedure computes the link integral of component ii and jj of writhe_pline
	using Banchoff's double dihedral sum formula. If ii and jj are the same, this 
        yields the writhe of this component. */

{
  int i, j;
  NUMPTR cur_writhe,ds,result,twopi,tmp,fourpi,tmp2;
  NUMPTR **np_buffer_ii;
  NUMPTR **np_buffer_jj;
  
  FILE   *outfile;
  double s_val, t_val, s_length, t_length, d_ds;
  int    point_cnt;

  /* As a prequel to the computation, establish a value for 2pi */

  tmp = (*double_to_np)(2.0);
  twopi = (*np_mul)(tmp,np_PI);  
  
  tmp2 = (*double_to_np)(4.0);
  fourpi = (*np_mul)(tmp2,np_PI);

  /* The first step is to convert the vertex buffer of writhe_pline
     into NUMVECTORS. First, we take care of the closed-polyline 
     issue. */

  int *edgebuf;
  edgebuf = calloc(writhe_pline->nc,sizeof(int));
  plc_edges(writhe_pline,edgebuf);

  np_buffer_ii = calloc(edgebuf[ii] + 1,sizeof(NUMPTR *));
  np_buffer_jj = calloc(edgebuf[jj] + 1,sizeof(NUMPTR *));

  for(i=0;i<edgebuf[ii]+1;i++) {

    np_buffer_ii[i] = vector_to_nv(writhe_pline->cp[ii].vt[i]);

  }

  for(i=0;i<edgebuf[jj]+1;i++) {

    np_buffer_jj[i] = vector_to_nv(writhe_pline->cp[jj].vt[i]);

  }

  /* Next, we do the dihedral sum. Note that if the polyline is closed,
     the last vertex is now a duplicate of the first, so we don't have to
     worry about missing a segment in the loop. 

     Also note that adjacent segments are coplanar, and hence can't
     contribute to the writhe. 

     If VERBOSE is on, we output a plot of the dihedral sum data to a 
     file, for later consumption by GNUPLOT. */

  cur_writhe = (*double_to_np)(0.0);

  if (VERBOSE) {

    fprintf(stderr,"np_writhe: Writing integrand data to np_integrand.dat...\n");
    outfile = fopen("np_integrand.dat","w");

    s_val = 0; t_val = 0; point_cnt = 0;
   
    for(i=1;i<edgebuf[ii]+1;i++) {

      s_length = plc_distance(writhe_pline->cp[ii].vt[i],writhe_pline->cp[jj].vt[i-1]);
      s_val += s_length;

      for(j=1;j<edgebuf[jj]+1;j++) {
	
	ds = np_dihedral_sum(np_buffer_ii[i],np_buffer_ii[i-1],       
			     np_buffer_jj[j],np_buffer_jj[j-1]);     
	cur_writhe = np_plusequal(cur_writhe,ds);
	
	t_length = plc_distance(writhe_pline->cp[jj].vt[j],writhe_pline->cp[jj].vt[j-1]);
	t_val   += t_length;
	
	d_ds = np_to_double(ds);
	d_ds *= 1.0/(4*3.141592653589793238462643383279); // 4 pi
	
	fprintf(outfile,"%g %g %g \n",t_val,s_val,d_ds);
	point_cnt++;

	np_integrand[i][j] = d_ds;
	
	(*kill_np)(ds);
	
      }
      
      t_val = 0;
      fprintf(outfile,"\n");
      
    }
    
    fclose(outfile);
    fprintf(stderr,"np_writhe: Output %d data points to file.\n\n",point_cnt);

    /* We now divide by 4 Pi (to normalize writhe by the area of the unit
       sphere). Notice that this differs from the ordinary case below: we
       have actually done twice as much work here by doing the full integral
       instead of only the lower half. */

     result = np_div(cur_writhe,twopi);
     result = np_div(result,tmp);   /* tmp still is equal to 2.0 */

  } else {

     for(i=3;i<edgebuf[ii]+1;i++) {

       for(j=1;j<edgebuf[jj]+1;j++) {
	
	 ds = np_dihedral_sum(np_buffer_ii[i],np_buffer_ii[i-1],       
			      np_buffer_jj[j],np_buffer_jj[j-1]);     
	 cur_writhe = np_plusequal(cur_writhe,ds);
	 (*kill_np)(ds);
	 
       }
      
     }

     /* We now divide by 4 Pi (to normalize writhe by the area of the unit
     sphere). */

     result = (*np_div)(cur_writhe,fourpi);

  }
   
  /* And free the spare variables */

  (*kill_np)(tmp); (*kill_np)(twopi); (*kill_np)(cur_writhe);
  (*kill_np)(tmp2); (*kill_np)(fourpi);
  
  for(i=0;i<edgebuf[ii]+1;i++) {

    for(j=0;j<3;j++) {

      (*kill_np)(np_buffer_ii[i][j]);

    }

    free(np_buffer_ii[i]);

  }

  free(np_buffer_ii);

  for(i=0;i<edgebuf[jj]+1;i++) {

    for(j=0;j<3;j++) {

      (*kill_np)(np_buffer_jj[i][j]);

    }

    free(np_buffer_jj[i]);

  }

  free(np_buffer_jj);

  /* Now we're really done. */

  return result;

}

double writhe(plCurve *writhe_pline,int precision,double *bycomponents)
     
     /* Procedure computes the writhe, in the current arithmetic model and
	precision, of writhe_pline, and returns the result as a double. 

	If bycomponents is not NULL, we also return the pairwise writhes 
	(or linking numbers) for each pair of components. We assume in this
	case that bycomponents is an array of writhe_pline->nc**2 doubles. 
	The writhe for pair (i,j) is stored in bycomponents[writhe_pline->nc*i + j];

	Why only a double? Well, usually, we are not trying to acheive more than
	16 digits of accuracy in the final result of a writhe computation; it
	is just that we must use more digits in the middle stages to reach the
	desired accuracy. */
     
{
  char   timestamp[100];
  double returned_writhe = 0;
  NUMPTR computed_writhe;

  double *writhegrid,ijwrithe;
  time_t start_comp, end_comp;

  if (bycomponents == NULL) {

    writhegrid = calloc(writhe_pline->nc * writhe_pline->nc,sizeof(double));

  } else {

    writhegrid = bycomponents;

  }

  (*init_arith_model)(precision);
  
  time(&start_comp);
  strftime(timestamp,sizeof(timestamp),"%a, %b %d (%I:%M %p)",localtime(&start_comp));
  
  if (VERBOSE == TRUE) {
    
    printf("Precision         : %d \n",  precision);
    printf("Start of run      : %s \n",  timestamp);
    
  }
  
  if (TIMING) {  /* We now quit */
    
    (*kill_arith_model)();
    exit(0);
    
  }
  
  // We are now ready to compute.

  int i,j;

  for(i=0;i<writhe_pline->nc;i++) {
    
    for(j=0;j<=i;j++) { 

      computed_writhe = np_writhe(writhe_pline,i,j);
      ijwrithe = (*np_to_double)(computed_writhe);

      writhegrid[writhe_pline->nc*i + j] = ijwrithe;
      writhegrid[writhe_pline->nc*j + i] = ijwrithe;

      if (i != j) { returned_writhe +=  2*ijwrithe; } else { returned_writhe += ijwrithe; }

      kill_np(computed_writhe);

    }

  }
  
  time(&end_comp);
  
  // We now report the results.
  
  if (VERBOSE == TRUE) {
    
    printf("writhe : %s \n\n",np_to_string(computed_writhe));
    strftime(timestamp,sizeof(timestamp),"%a, %b %d (%I:%M %p)",localtime(&end_comp));
    
    printf("End of run        : %s\n",timestamp);
    printf("Total time        : %d sec.\n\n",(int)(difftime(end_comp,start_comp)));
    
  }
  
  /* Now we shut down the model. */
  
  kill_arith_model();

  if (bycomponents == NULL) { free(writhegrid); }
  
  return returned_writhe;
  
}

int arithmodel_ok()

     /* Procedure performs a series of confidence tests on the current 
	arithmetic model, all of which are designed to make sure that the
	computation engine is producing reasonable results. 

        If VERBOSE is turned on, produces a report of tests passed and 
        failed as we proceed. Otherwise, simply return TRUE or FALSE,
        depending on whether the model passes its confidence tests. */

{
  char *PIstring,*arithstring[3];
  NUMPTR one;
  char   printstring[256];

  char pistring[128] = {"3.141592653589793238462643383279502884197169399375"\
		   "1058209749445923078164"}; /* 70 digits of PI */

  char sqrt2[128] = {"1.41421356237309504880168872420969807856967187537694"\
		"80731766797379907324"};     /* 70 digits of sqrt 2 */

  double x;
  NUMPTR X,Y,Z,W;
  NUMPTR tmp[5],zero;
  int    i;

  NUMVECTOR nvA,nvB;
  NUMPTR    *nvRes;

  init_arith_model(DEC_PRECISION);

  /* Test constants. */

  PIstring = (*np_to_string)(np_PI);
  
  if (strncmp(pistring,PIstring,DEC_PRECISION+1) != 0) { /* Strings don't match. */

    if (VERBOSE) {

      sprintf(printstring,"arithmodel_ok: Model constant PI != correct value.\n"\
	      "               PI       = (%%.%ds) \n"\
	      "               True val = (%%.%ds) \n\n",DEC_PRECISION,DEC_PRECISION);

      fprintf(stderr,printstring,PIstring,pistring);

    }
      
    return FALSE;

  } else if (VERBOSE) {

    sprintf(printstring,"arithmodel_ok: Model constant PI == correct value.\n"\
	    "               PI       = (%%.%ds) \n"\
	    "               True val = (%%.%ds) \n\n",DEC_PRECISION,DEC_PRECISION);

    fprintf(stderr,printstring,PIstring,pistring);
    fprintf(stderr,"arithmodel_ok: 1st test PASSED.\n");

  }

  /* Now test np_MINUS to make sure that it's accurate. */

  x = (*np_to_double)(np_MINUS);

  if (fabs(x + 1) > 1e-14) {

    if (VERBOSE) {

      fprintf(stderr,"arithmodel_ok: Model constant -1 != correct value.\n"\
	      "               np_MINUS = %16g.\n",x);
      
    }

    return FALSE;

  } else if (VERBOSE) {

    fprintf(stderr,"arithmodel_ok: Model constant -1 == correct value.\n"\
	    "               np_MINUS = %16g.\n\n",x);
    fprintf(stderr,"arithmodel_ok: 2nd test PASSED.\n");

  }

  /* Now we test double_to_np and np_to_double. */

  x = 1.0/7.0;
  X = (*double_to_np)(x);
  x = (*np_to_double)(X);

  if (fabs(x - 1.0/7.0) > 1e-14) {

    if (VERBOSE) {

      sprintf(printstring,"arithmodel_ok: np_to_double or double_to_np have failed.\n"\
	                  "               original value      : %%15g \n"\
	                  "               converted to NUMPTR : %%s\n"\
	      "               converted to double : %%15g \n");

      fprintf(stderr,printstring,1.0/7.0,np_to_string(X),x);

    }

    return FALSE;

  } else if (VERBOSE) {

      sprintf(printstring,"arithmodel_ok: np_to_double and double_to_np test ok.\n"\
	                  "               original value      : %%15g \n"\
	                  "               converted to NUMPTR : %%s\n"\
	      "               converted to double : %%15g \n");

      fprintf(stderr,printstring,1.0/7.0,(*np_to_string)(X),x);
  
      fprintf(stderr,"\narithmodel_ok: 3rd test PASSED.\n");

  } 

  (*kill_np)(X);

  /* Now we test the elementary arithmetic functions. */

  X = (*double_to_np)(-3.0);
  Y = (*double_to_np)(-55.67);
  Z = (*double_to_np)(-10000.0);
  W = (*double_to_np)(-0.0001);

  tmp[0] = (*np_add)(X,Y);
  tmp[1] = (*np_mul)(tmp[0],Z);
  tmp[2] = (*np_div)(tmp[0],W);

  arithstring[0] = (*np_to_string)(tmp[1]);
  arithstring[1] = (*np_to_string)(tmp[2]);

  if (strncmp(arithstring[0],arithstring[1],DEC_PRECISION-1) != 0) {  /* The strings _differ_ */

    if (VERBOSE) {

      sprintf(printstring,"arithmodel_ok: add, mul, or div failed.\n"\
	                  "               mul chain : %%s \n"\
	                  "               div chain : %%s \n"\
	                  "               expected  : 586700 \n");

      fprintf(stderr,printstring,arithstring[0],arithstring[1]);
    
    }

    return FALSE;

  } else if (VERBOSE) {

    sprintf(printstring,"arithmodel_ok: add, mul, and div passed.\n"\
	                  "               mul chain : %%s \n"\
	                  "               div chain : %%s \n"\
	                  "               expected  : 586700 \n");

    fprintf(stderr,printstring,arithstring[0],arithstring[1]);

    fprintf(stderr,"\narithmodel_ok: 4th test PASSED.\n");

  }

  /* Now running many arithmetic operations to test memory usage. */

  one = double_to_np(1.0);
  zero = double_to_np(0.0);

  for (i=0;i<10000;i++) {

    tmp[1]= np_mul(one,one);
    kill_np(tmp[1]);
    
  }

  for (i=0;i<10000;i++) {

    tmp[1] = np_div(one,one);
    kill_np(tmp[1]);

  }

  for (i=0;i<10000;i++) {

    tmp[1] = np_add(one,zero);
    kill_np(tmp[1]);

  }

  for (i=0;i<10000;i++) {

    tmp[1] = double_to_np(1.0);
    kill_np(tmp[1]);

  }

  for (i=0;i<10000;i++) {

    tmp[1] = np_acos(zero);
    kill_np(tmp[1]);

  }

  for (i=0;i<10000;i++) {

    tmp[1] = np_sqrt(one);
    kill_np(tmp[1]);

  }
  
  for (i=0;i<10000;i++) {

    np_leq(zero,one);

  }

  tmp[1] = one;

  for (i=0;i<10000;i++) {

    tmp[1] = np_minusequal(tmp[1],zero);

  }

  for(i=0;i <10000; i++) {

    tmp[1] = np_plusequal(tmp[1],zero);

  }

  nvA[0] = double_to_np(1.0); 
  nvA[1] = double_to_np(2.0);
  nvA[2] = double_to_np(3.0);

  nvB[0] = double_to_np(1.0);
  nvB[1] = double_to_np(3.0);
  nvB[2] = double_to_np(4.0);

  for(i=0; i < 10000;i++) {
    
    nvRes = nv_diff(nvA,nvB);
    kill_nv(nvRes);
    free(nvRes);

  }

  for(i=0; i < 10000;i++) {
    
    nvRes = nv_cross(nvA,nvB);
    kill_nv(nvRes);
    free(nvRes);

  }

  for(i=0; i < 10000;i++) {
    
    tmp[0] = nv_dot(nvA,nvB);
    kill_np(tmp[0]);

  }

  for(i=0; i < 10000;i++) {

    nv_normalize(nvA);

  }


  (*kill_np)(X); 
  (*kill_np)(Y); 
  (*kill_np)(Z); 
  (*kill_np)(W); 
  
  (*kill_np)(tmp[0]);
  (*kill_np)(tmp[1]);
  (*kill_np)(tmp[2]);

  free(arithstring[0]);
  free(arithstring[1]);

  /* We now test the acos function by writing pi = 3 acos (1/2). */

  X = (*double_to_np)(3.0);
  Y = (*double_to_np)(1.0);
  Z = (*double_to_np)(2.0);

  W = (*np_div)(Y,Z);

  tmp[0] = (*np_acos)(W);
  tmp[1] = (*np_mul)(X,tmp[0]);

  arithstring[0] = (*np_to_string)(tmp[1]);

  if (strncmp(pistring,arithstring[0],DEC_PRECISION+1) != 0) { 
  
    if (VERBOSE) {

      sprintf(printstring,"arithmodel_ok: acos computation of Pi failed.\n"\
	                  "               acos pi     = (%%.%ds) \n"\
	      "               True val = (%%.%ds) \n\n",DEC_PRECISION,DEC_PRECISION);

      fprintf(stderr,printstring,arithstring[0],pistring);

    }
      
    return FALSE;

  } else if (VERBOSE) {

    sprintf(printstring,"arithmodel_ok: acos computation of pi passed.\n"\
	    "               acos pi  = (%%.%ds) \n"\
	    "               True val = (%%.%ds) \n\n",DEC_PRECISION,DEC_PRECISION);

    fprintf(stderr,printstring,arithstring[0],pistring);
    fprintf(stderr,"arithmodel_ok: 5th test PASSED.\n");

  }
  
  (*kill_np)(X); (*kill_np)(Y); (*kill_np)(Z); (*kill_np)(W);
  free(arithstring[0]);

  /* 6th test: We make sure sqrt is working by computing sqrt(2). */

  X = (*double_to_np)(2.0);
  Y = (*np_sqrt)(X);
  
  arithstring[0] = (*np_to_string)(Y);

  if (strncmp(sqrt2,arithstring[0],DEC_PRECISION+1) != 0) { /* Strings don't match. */

    if (VERBOSE) {

      sprintf(printstring,"arithmodel_ok: sqrt(2) computation failed.\n"\
	      "               sqrt(2)     = (%%.%ds) \n"\
	      "               True val = (%%.%ds) \n\n",DEC_PRECISION,DEC_PRECISION);

      fprintf(stderr,printstring,arithstring[0],sqrt2);

    }
      
    return FALSE;

  } else if (VERBOSE) {

    sprintf(printstring,"arithmodel_ok: sqrt(2) computation passed.\n"\
	    "               sqrt(2)  = (%%.%ds) \n"\
	    "               True val = (%%.%ds) \n\n",DEC_PRECISION,DEC_PRECISION);

    fprintf(stderr,printstring,arithstring[0],sqrt2);
    fprintf(stderr,"arithmodel_ok: 6th test PASSED.\n");

  }
  
  (*kill_np)(X); (*kill_np)(Y);
  free(arithstring[0]);
  
  /* 7th test. Testing np_leq. */

  X = (*double_to_np)(1.0);
  Y = (*double_to_np)(1.0001);
  Z = (*double_to_np)(1.0);

  arithstring[0] = (*np_to_string)(X);
  arithstring[1] = (*np_to_string)(Y);
  arithstring[2] = (*np_to_string)(Z);

  if ((*np_leq)(Y,X) || !(*np_leq)(X,Z) || !(*np_leq)(Z,X) || !(*np_leq)(X,Y)) {

    if (VERBOSE) {

      fprintf(stderr,"arithmodel_ok: one of np_leq tests failed. \n"\
	      "               %s <= %s : %d (expected 0) \n"\
	      "               %s <= %s : %d (expected 1) \n"\
	      "               %s <= %s : %d (expected 1) \n"\
	      "               %s <= %s : %d (expected 1) \n",
	      arithstring[1],arithstring[0],(*np_leq)(Y,X),
	      arithstring[0],arithstring[2],(*np_leq)(X,Z),
	      arithstring[2],arithstring[0],(*np_leq)(Z,X),
	      arithstring[0],arithstring[1],(*np_leq)(X,Y));
      
    }

    return FALSE;

  } else if (VERBOSE) {

    fprintf(stderr,"arithmodel_ok: np_leq tests passed. \n"\
	    "               %s <= %s : %d (expected 0) \n"\
	    "               %s <= %s : %d (expected 1) \n"\
	    "               %s <= %s : %d (expected 1) \n"\
	    "               %s <= %s : %d (expected 1) \n",
	    arithstring[1],arithstring[0],(*np_leq)(Y,X),
	    arithstring[0],arithstring[2],(*np_leq)(X,Z),
	    arithstring[2],arithstring[0],(*np_leq)(Z,X),
	    arithstring[0],arithstring[1],(*np_leq)(X,Y));
    
    fprintf(stderr,"arithmodel_ok: 7th test PASSED.\n");

  }
  
  kill_arith_model();

  /* Tests are now complete. */

  if (VERBOSE) {

    fprintf(stderr,"\narithmodel_ok: All model tests PASSED.\n\n");
    
  }

  return TRUE;

}
  
	    
	    

    

