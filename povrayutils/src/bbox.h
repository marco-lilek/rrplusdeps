/*

  bbox.h : Prototypes for procedures to read and write .bbox files. 

*/

#ifndef BBOX_H__
#define BBOX_H__ 1

#include "config.h"

#ifdef HAVE_STDBOOL_H
#include<stdbool.h>
#endif

#include "plCurve.h"

void write_bbox(plc_vector l,plc_vector u,FILE *outfile);
bool read_bbox(plc_vector *l,plc_vector *u,FILE *infile);

#endif
