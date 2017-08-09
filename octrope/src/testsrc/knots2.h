#include<iostream>
#include<assert.h>
#include"point.h"
#include"vector.h"
#include"vectors.h"

#ifndef KNOTS_H
#define KNOTS_H

using namespace std;

class knots{
public:
  knots();
  knots(int numberOfSides);
  /* knots(int p, int q, double bigRadius,
	double littleRadius, int numberOfVertices); */
  knots(const knots & source);
  ~knots();
  
  const point & operator [](const int & index) const;
  point & operator [](const int & index);
  
  knots & operator =(const knots & source); 

  vectors lenvtr;
  vector side;
  vectors angle;

  void Compute();
  void rotatexy(const double & theta);
  void rotateyz(const double & angle);
  void rotaterawyz(const double & cosangle);
  void rotaterawnegyz(const double & cosangle);
  void rotatezx(const double & phi);
  void rotateyz(const double & rho, const int & i, const int & j);
  void BeShaken(const vector & shaker);

  void Display();
  double injrad() const;
  double Length() const;
  double RopeLength() const;
  double Thickness() const;
  double MinRad() const;
  double InteriorInterior() const;
  double VertexInterior() const;
  double VertexVertex() const;
  void SideData() const;
  void AngleData() const;
  void CurvatureData() const;
  void CurvatureSpreadsheet() const;
  double minRadFactor() const;
  double ChangeMinRadFactor(double newfactor);
  void CriticalRelations(); 

  int setSize(int numberOfSides);
  int vnum() const;
  double Epsilon() const;
  int find_rate(const int & i, const int & j, 
		double & a, double & b, double & d) const;
  int IsDC(const int & i, const int & j, const double & a) const;

  friend istream & operator >> (istream & in, knots & source);

private:
  double R[4];
  double minRad_factor;
  double injRadius;
  double epsilon;
  double polyLength;
  void ThicknessData();
  vector vertex;
  int numSides;
};

ostream & operator << (ostream & out, const knots & source);
int operator ==(const knots & a, const knots & b);
knots operator + (const knots & knot1, const knots & knot2);
knots operator * (double num, const knots & source);
int operator !=(const knots & a, const knots & b);

#endif
