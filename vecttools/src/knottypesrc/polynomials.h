/*****

	polynomials.h : Definitions for the polynomials source file.


********/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_STDBOOL_H
#include <stdbool.h>
#endif

#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif

#include"../utilib/ordie.h"

typedef struct monostruct {

	double coeff;
	int l;
	int m;

} monomial;

monomial *lmpoly_to_polynomial(char *lmoutput,int *nmonomials);
char     *lmpoly_to_mathematica(char *lmpoly_output);
char     *lmpoly_to_latex(char *lmpoly_output);

char *poly_to_mathematica(monomial *poly, int *nmonos);
char *poly_to_latex(monomial *poly, int *nmonos);

monomial *product_polynomial(monomial *pA, int nA, monomial *pB, int nB, int *nProduct);
char     *polynomial_to_lmpoly(monomial *poly,int nmonos);

char *lmpoly_check(char *lmpoly_output);

