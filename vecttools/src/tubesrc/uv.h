/*

   uv.h : Procedures for appending u-v texture space coordinates to a plSurf. 

*/

#ifndef __UV_H__
#define __UV_H__ 1

struct uvbuf {

  int verts;
  plc_vector *uv;  /* uv[i].c[0] = u, uv[i].c[1] = v, uv[i].c[2] = 0 */ 

};

void write_uvfile(struct uvbuf *uvb,FILE *outfile);
int  load_uvfile(struct uvbuf *uvb,FILE *infile);

#endif
