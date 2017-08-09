#include<iostream>
#include<assert.h>
#include<math.h>

#ifndef POINT_H
#define POINT_H

using namespace std;

class point{
public:
  point();
  point(double a, double b, double c);
  point(const point & source);

  double & operator [] (int index);
  const double & operator [] (int index) const;

  point & operator = (const point & source);

  double norm();

  double coord[3];
};

ostream & operator << (ostream & out, const point & source);
point operator +(point a, point b);
point operator -(point a, point b);
point operator -(point a);
point operator *(double x, point b);
int operator !=(point a, point b);
int operator ==(point a, point b);
double operator *(point a, point b);
point operator &(point a, point b);

#endif
