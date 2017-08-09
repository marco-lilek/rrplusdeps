#include<iostream>
#include<assert.h>
#include"point.h"
#include"vector.h"
#include"knots2.h"
// #include"tcknots.h"

#ifndef SHORTKNOT_H
#define SHORTKNOT_H

using namespace std;

class shortknot{
public:
  shortknot();
  shortknot(unsigned int numberOfSides);
  /* shortknot(unsigned int p, unsigned int q, double bigRadius,
	double littleRadius, int numberOfVertices); */
  shortknot(unsigned int numberOfSides, double fillvalue);
  shortknot(const shortknot & source);
  shortknot(const knots & source);
  /*  shortknot(const tcknots & source); */
  ~shortknot();
  
  const point & operator [](int index) const;
  point & operator [](int index);

  unsigned int setSize(unsigned int numberOfSides);
  unsigned int vnum() const;
  void clear();

  friend istream & operator >> (istream & in, shortknot & source);

private:
  vector vertex;
  unsigned int numSides;
};

ostream & operator << (ostream & out, const shortknot & source);
unsigned int operator ==(shortknot a, shortknot b);
shortknot operator + (const shortknot & knot1, const shortknot & knot2);
shortknot operator * (double num, const shortknot & source);
unsigned int operator !=(shortknot a, shortknot b);

#endif
