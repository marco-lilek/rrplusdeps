#include<iostream>
#include<assert.h>
#include"point.h"

#ifndef VECTOR_H
#define VECTOR_H

using namespace std;

class vector{
public:
  vector();
  vector(int numElts);
  vector(const vector & source);
  ~vector();

  point & operator[] (int index);
  const point & operator [] (int index) const;

  vector & operator =(const vector & source);

  int length () const;

  int setSize(int numberOfElements);

protected:
  point * data;
  int size;
};

vector operator *(const double & scalar, const vector & source);

#endif
