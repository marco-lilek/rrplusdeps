/* This program minimizes polygonal ropelength in two ways:
   
   1) keeping the side length fixed (+- computational error)
   2) randomly shaking vertices with no regard to edge lengths

   Each run is a technique of handling the data determined by such shaking.

   Compiling Instructions:
     CC -c *.C
     CC toros.cc *.o -lm
   
   Eric Rawdon 2/98
*/


#include<iostream>
#include<math.h>
#include<fstream>
#include<string>
#include<assert.h>
#include"knots2.h"
#include"vector.h"
#include"point.h"
#include"vectors.h"
// #include"CPstring.h"
#include"shortknot.h"
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>
#include"miniball.H"
#include<gsl/gsl_math.h>
#include<gsl/gsl_eigen.h>

/* Function Prototypes */
void Shake(knots & knot, double magnitude);
void ACWN(const knots & knot, double & acn, double & awn);
double ComputeMaxAngle(const knots & knot, const int & i, const int & j);
void Run(knots & knot, int runs);
void ConjugateRun(knots & knot, int runs);
void PersistentRun(knots & knot, int runs);
void AnnealingRun(knots & knot, double temp, 
		  const double & dissipation, int runs);
void VertexRun(knots & knot, int runs);
void EqRun(knots & knot, int runs);
void Tidy(knots & knot, int runs);
void EqConjugateRun(knots & knot, int runs);
void EqPersistentRun(knots & knot, int runs);
void EqAnnealingRun(knots & knot, double temp, 
		  const double & dissipation, int runs);
void DoubleSides(knots & knot);
void Menu();
void Load(knots & knot, string & file);
void LoadMillett(knots & knot, const string & file);
void Save(const knots & knot, const string & file);
void SaveMillett(const knots & knot, const string & file);
void ShowVertexInfo(const knots & knot);
knots MakeEquilateral(const knots & knot);
knots MakeEquilateral2(knots & knot);
void Normalize2(knots & knot, const int & i, const int & j);
void PointsNormalize(knots & knot, const point & v1, const point & v2);
void Resize(knots & knot, const double & newsize);
void Standardize(knots & knot, double & d1, double & d2, double & d3);
void Standardize2(knots & knot, double & d1, double & d2, double & d3);
void Standardize3(knots & knot, double & d1, double & d2, double & d3);
void AcnAnneal(knots & knot, const double & temperature, int runs);
void HowThickSmooth(const knots & knot);
void HowThickSmooth2(const knots & knot);
void HowThickSmoothFactor(const knots & knot);
void IsDCincharge(const knots & knot);
void CheckCriticality(const knots & knot);
void CheckCriticality2(const knots & knot);
int NewFindRate(const knots & knot, const int & i, const int & j, double & d);
double min(const knots & knot, const int & i, const point & p, double & ratio);
int IsDC(const knots & knot, const int & i, const int & j, const double & a);
void Torsion(knots & knot);
void Torsion2(knots & knot);
void ConvexHull(knots & knot);
void QHull(knots & knot);
void Hull2(knots & knot, double & volume, double & surfacearea);
void EqCriticalRun(knots & knot, const int & runs, const double & anglechange);
void AutomateFindCrit(knots & knot);
void NewBoxDimensions(knots & knot);
void SkinnyMinBoxDimensions(knots & knot, double &, double &, double &);
void SkinnyMaxBoxDimensions(knots & knot, double &, double &, double &);
void IsItCrit(knots & knot);
void CritFind(knots & knot);
void NonEqCritFind(knots & knot);
void Distance(const knots & knot);
void IsItEq(knots & knot);
void DefaultPosition(knots & knot);
void Mobius(knots & knot);
void Report(knots & knot, const string & file);
void CurvaturePlot(knots & knot);
void RadiusOfGyration(const knots & knot);
int  MiniBall(const knots & knot);
void MakeThicknessOne(knots & knot);
double Asphericity(const knots & knot);
double piatekasphprol(const knots & knot);
double MDEnergy(const knots & knot);
double SymmetryScore(const knots & knot);

int main()
{
  knots knot;

  string file;
 
  /* Load an initial file */
  cout << "File name: ";
  cin >> file;
  Load(knot,file);
  
  int runs;
  double temp, dissipation;
  string command;

  do{
    cout << "Action: ";
    cin >> command;
    
    if(command == "report")
      {
	Report(knot,file);
      }

    else if(command == "radgyr")
      {
	RadiusOfGyration(knot);
      }

    else if(command == "r")
      {
	cout << "# runs: ";
	if(cin >> runs)
	  Run(knot,runs);
      }

    else if(command == "er")
      {
	cout << "# runs: ";
	if(cin >> runs)
	  EqRun(knot,runs);
      }
    
    else if(command == "a")
      {
	cout << "temperature dissipation runs: ";
	if( (cin >> temp) && (cin >> dissipation) && (cin >> runs) )
	  AnnealingRun(knot, temp, dissipation, runs);
      }

    else if(command == "ea")
      {
	cout << "temperature dissipation runs: ";
	if( (cin >> temp) && (cin >> dissipation) && (cin >> runs) )
	  EqAnnealingRun(knot, temp, dissipation, runs);
      }
    
    else if(command == "p")
      {
	cout << "# runs: ";
	if(cin>>runs)
	  PersistentRun(knot,runs);
      }
    
    else if(command == "ep")
      {
	cout << "# runs: ";
	if(cin>>runs)
	  EqPersistentRun(knot,runs);
      }
    
    else if(command == "c")
      {
	cout << "# runs: ";
	if(cin>>runs)
	  ConjugateRun(knot,runs);
      }
    
    else if(command == "ec")
      {
	cout << "# runs: ";
	if(cin>>runs)
	  EqConjugateRun(knot,runs);
      }

    else if(command == "mobius")
      {
	Mobius(knot);
      }

    else if(command == "findcrit")
      {
	AutomateFindCrit(knot);
      }

    else if(command == "eqcrit")
      {
	cout << "Number of iterations: ";
	if(cin >> runs){
	  cout << "Amount to perturb the angle: ";
	  double anglechange;
	  if(cin >> anglechange)
	    EqCriticalRun(knot,runs,anglechange);
	}
      }

    else if(command == "acnrun")
      {
	cout << "Temperature: ";
	if(cin>>temp)
	  {
	    cout << "# runs: ";
	    if(cin>>runs)
	      AcnAnneal(knot,temp,runs);
	  }
      }

    else if(command == "scale")
      {
	cout << "New length: ";
	if(cin>>temp)
	  {
	    Resize(knot,temp);
	  }
      }
  
    else if(command == "mt1")
      {
	MakeThicknessOne(knot);
      }
  
    else if(command == "asphericity")
      {
	cout << "Asphericity: " << Asphericity(knot) << endl;
      }
  
    else if(command == "symmetry"){
      SymmetryScore(knot);
    }

    else if(command == "energy")
      {
	cout << "Normalized MD-Energy: " << MDEnergy(knot) << endl;
      }
  
    else if(command == "tidy")
      {
	cout << "# runs: ";
	if(cin>>runs)
	  Tidy(knot,runs);
      }
    
    else if(command == "v")
      {
	cout << "# runs: ";
	if(cin>>runs)
	  VertexRun(knot,runs);
      }

    else if(command == "acwn")
      {
	double acn, awn, totalcurvature = 0;
	ACWN(knot, acn, awn);
	cout.setf(ios::fixed);
	cout.precision(6);
	cout << "Average crossing number is " << acn << endl;
	cout << "Average writhe is " << awn << endl;
	int top = knot.vnum();
	for(int i=0; i<top; i++)
	  totalcurvature += knot.angle[i];
	cout << "Total Curvature is " << totalcurvature << endl;
      }

    else if(command == "torsion")
      {
	Torsion(knot);
      }

    else if(command == "torsion2")
      {
	Torsion2(knot);
      }

    else if(command == "factor")
      {
	double factor;
	cout << "New factor must lie between 0 and 1.  Default is 1.\n";
	cout << "New factor: ";
	if(cin >> factor)
	  knot.ChangeMinRadFactor(factor);
      }

    else if(command == "miniball")
      MiniBall(knot);

    else if(command == "smooth")
      HowThickSmooth(knot);

    else if(command == "smooth2")
      HowThickSmooth2(knot);

    else if(command == "isdc")
      IsDCincharge(knot);

    else if(command == "check")
      CheckCriticality(knot);

    else if(command == "check2")
      CheckCriticality2(knot);

    else if(command == "q" || command == "Q")
      {
	break;
      }

    else if(command == "dt")
      knot.Display();

    else if(command == "ds")
      knot.SideData();

    else if(command == "da")
      knot.AngleData();

    else if(command == "dc")
      knot.CurvatureData();

    else if(command == "cp")
      CurvaturePlot(knot);

    else if(command == "sdc")
      knot.CurvatureSpreadsheet();

    else if(command == "dv")
      cout << knot << endl;

    else if(command == "vertex")
      ShowVertexInfo(knot);

    else if(command == "load")
      {
	cout << "File name: ";
	if(cin >> file)
	  Load(knot,file);
      }

    else if(command == "loadm")
      {
	cout << "File name: ";
	if(cin >> file)
	  LoadMillett(knot,file);
      }	
    
    else if(command == "ls")
      {
	system("ls");
      }

    else if(command == "system")
      {
	string help;
	char comma;
	cout << "Type system command: ";
	cin.get(comma);
	getline(cin,help);
	system(help.c_str());
      }

    else if(command == "save")
      {
	cout << "File name: ";
	if(cin >> file)
	  Save(knot,file);
      }

    else if(command == "savem")
      {
	cout << "File name: ";
	if(cin >> file)
	  SaveMillett(knot,file);
      }

    else if(command == "double")
      {
	DoubleSides(knot);
	knot.Compute();
	cout << "Number of sides is now " << knot.vnum() << endl;
      }

    else if(command == "box")
      {
	double d1,d2,d3;
	Standardize(knot,d1,d2,d3);
      }

    else if(command == "minbox")
      {
	double length, width, height;
	SkinnyMinBoxDimensions(knot,length,width,height);
      }

    else if(command == "maxbox")
      {
	double length, width, height;
	SkinnyMaxBoxDimensions(knot,length,width,height);
      }

    else if(command == "ellipse")
      {
	double d1,d2,d3;
	Standardize2(knot,d1,d2,d3);
      }

    else if(command == "equilize")
      {
	MakeEquilateral2(knot);
      }
    
    else if(command == "convex")
      {
	ConvexHull(knot);
      }

    else if(command == "qhull")
      {
	QHull(knot);
      }

    else if(command == "hull")
      {
	double volume, surfacearea;
	Hull2(knot,volume,surfacearea);
      }

    else if(command == "newbox")
      {
	NewBoxDimensions(knot);
      }

    else if (command == "isitcrit")
      {
	IsItCrit(knot);
      }

    else if (command == "findrealcrit")
      {
	CritFind(knot);
      }

    else if (command == "noneqcrit")
      {
	NonEqCritFind(knot);
      }

    else if (command == "distance")
      {
	Distance(knot);
      }

    else if (command == "isiteq")
      {
	IsItEq(knot);
      }

    else if (command == "seedc")
      {
	knot.CriticalRelations();
      }

    else if (command == "default")
      {
	DefaultPosition(knot);
      }

    else
      {
	Menu();
      }
  }while(1);

  Save(knot, "boogie");
  return 1;
}      
      

void Menu()
{
  cout << "Menu\n";
  cout << "-------------------------------\n";
  cout << "Run: r runs\n";
  cout << "Equilateral Run: er runs\n";
  cout << "Anneal: a temp dissipation runs\n";
  cout << "Equilateral Anneal: ea temp dissipation runs\n";
  cout << "Conjugate: c runs\n";
  cout << "Equilateral Conjugate: ec runs\n";
  cout << "Conjugate Vertex: v runs\n";
  cout << "Persistent: p runs\n";
  cout << "Equilateral Persistent: ep runs\n";
  cout << "Double Number of Sides: double\n";
  cout << "Compute ACN and AWN: acwn\n";
  cout << "Equilateral ACN anneal: acnrun temperature runs\n";
  cout << "Equilateral Critical Run: eqcrit runs change\n";
  cout << "Find nearby critical conformation: findcrit\n";
  cout << "Change minRad factor: factor\n";
  cout << "Make Equilateral: equilize\n";
  cout << "Display Thickness Information: dt\n";
  cout << "Display Vertices: dv\n";
  cout << "Display Side Information: ds\n";
  cout << "Display Angle Information: da\n";
  cout << "Display Curvature Information: dc\n";
  cout << "Display Vertex Information: vertex\n";
  cout << "Display Bounds on an Inscribed Smooth Knot: smooth\n";
  cout << "Display Box Dimensions: box\n";
  cout << "Display Ellipsoid Dimensions: ellipse\n";
  cout << "Display torsion: torsion\n";
  cout << "Check the criticality numbers (w/ tc): check\n";
  cout << "Check the criticality numbers (no tc): check2\n";
  cout << "See doubly critical data: seedc\n";
  cout << "Do a Mobius transformation to the knot: mobius\n";
  cout << "Load knot: save filename\n";
  cout << "Save knot: save filename\n";
  cout << "Save knot Mathematica style: savem filename\n";
  cout << "quit: q\n\n";
}

/* --------------------------------------------------------------------- */
/*            The runs below do not conserve edge lengths                */
/* --------------------------------------------------------------------- */


/* Random perturbation of the vertices with a maximum coordinate 
   change = magnitude */
void Shake(knots & knot, double magnitude)
{
  int top = knot.vnum(), i, j;
  magnitude *= 2;

  for(i=0; i < top; i++)
    for(j=0; j<3; j++)
      knot[i][j] = knot[i][j] + magnitude * (drand48()-.5);
}

/* Normal run with random perturbations */
void Run(knots & knot, int runs)
{
  knots tempknot;

  for(; runs > 0; runs--)
    {
      tempknot = knot;
      Shake(tempknot, knot.Epsilon() * .1);
      tempknot.Compute();
      /* compareThickness(knot, tempknot); */
      if (tempknot.RopeLength() <= knot.RopeLength() )
	knot = tempknot;
    }

  cout << "Final knot data:" << endl;
  knot.Display();
}

/* Conjugate run attempts successful shakes until unsuccessful */
void ConjugateRun(knots & knot, int runs)
{
  int i, j;
  int top = knot.vnum();
  vector shaker(top);

  for(; runs > 0; runs--)
    {
      knots tempknot;
      double magnitude = knot.Epsilon() * .1;

      /* fill shaker with a random knot in case shake is successful*/
      for(i=0; i < top; i++)
	for(j=0; j<3; j++)
	  shaker[i][j] = magnitude * (drand48()-.5);

      /* compare the ropelengths and keep trying the same vertices as long
	 as it keeps working */
      do{
	tempknot = knot;
	knot.BeShaken(shaker);
	knot.Compute();
      } while (knot.RopeLength() <= tempknot.RopeLength());

      knot = tempknot;
    }

  cout << "Final knot data:" << endl;
  knot.Display();
}

/* Persistant attempts smaller perturbations in the shake direction
   until it gives up or achieves a successful ropelength reduction */
void PersistentRun(knots & knot, int runs) 
{ 
  int i,j;
  int top = knot.vnum();
  vector shaker(top);
  
  for(; runs > 0; runs--)
    {
      knots tempknot(knot);
      double magnitude = knot.Epsilon() * .1;
      /* fill shake with a random knot */
      for(i=0; i < top; i++)
	for(j=0; j<3; j++)
	  shaker[i][j] = magnitude * (drand48()-.5);

      knot.BeShaken(shaker);
      knot.Compute();

      /* compare the ropelengths, if it doesn't work, try a smaller amount */
      if(knot.RopeLength() <= tempknot.RopeLength())
	{ 
	  do{
	    tempknot = knot;
	    knot.BeShaken(shaker);
	    knot.Compute();
	  } while (knot.RopeLength() <= tempknot.RopeLength());
	}
      else
	{
	  int h=0;
	  do{
	    knot = tempknot;
	    h++;
	    shaker = .5*shaker;
	    knot.BeShaken(shaker);
	    knot.Compute();
	  } while( (knot.RopeLength() > tempknot.RopeLength()) && (h<11) );
	}	  

      knot = tempknot;
    }

  cout << "Final knot data:" << endl;
  knot.Display();
}

/* Annealing always keeps successful shakes and sometimes keeps shakes
   which increase the ropelength */
void AnnealingRun(knots & knot, double temp, 
		  const double & dissipation, int runs)
{
  knots shortest(knot);
  double shortestLength = knot.RopeLength();
  double qtemp = 3000/temp;   /* try to speed things up */

  for(; runs > 0; runs--)
    {
      knots tempknot(knot);
      double tempRopeLength = knot.RopeLength();

      Shake(knot, knot.Epsilon() * .1);
      knot.Compute();
      /* compareThickness(knot, tempknot); */
      if (knot.RopeLength() > tempRopeLength)
	{
	  double probability = drand48();
	  /* percent is percent change from old rope length and is negative */
	  double percent = (tempRopeLength-knot.RopeLength())/tempRopeLength;
	  if(exp(percent*qtemp)<drand48())
	    {
	      knot = tempknot;
	    }
	}
      /* shortest keeps the knot with the lowest ropelength achieved */
      else if(knot.RopeLength() <= shortestLength)
	{
	  shortestLength = knot.RopeLength();
	  shortest = knot;
	}      
    }

  string filename;

  cout << "Save the shortest to file: ";
  cin >> filename;

  Save(shortest, filename);
  cout << "Final knot data:" << endl;
  cout << "Shortest was " << shortest.RopeLength() << endl;
  knot.Display();
}

/* Vertex tries one vertex at a time in a conjugate fashion */
void VertexRun(knots & knot, int runs)
{
  /* this is actually a conjugate vertex run */
  int i, j;
  int top = knot.vnum();
  point shaker;

  for(; runs > 0; runs--)
    {
      knots tempknot;
      for(i=0; i < top; i++)
	{
	  double magnitude = knot.Epsilon() * .1;

	  for(j=0; j<3; j++)
	    shaker[j] = magnitude * (drand48()-.5);
	  
	  /* compare the ropelengths and keep trying the same vertices as long
	     as it keeps working */
	  do{
	    tempknot = knot;
	    knot[i] = knot[i] + shaker;
	    knot.Compute();
	  } while (knot.RopeLength() <= tempknot.RopeLength());
	  
	  knot = tempknot;
	}
    }
  
  cout << "Final knot data:" << endl;
  knot.Display();
}

/* ----------------------------------------------------------------------- */
/* Equilateral dividing line, all functions below conserve side lengths +- */
/* ----------------------------------------------------------------------- */


/* translates knot so that vertex_i is at the origin and vertex_j lies
   on the positive x-axis (y and z coords are 0)*/
void Normalize(knots & knot, const int & i, const int & j)
{
  int k, top = knot.vnum();

  /* Translate */
  point newOrigin = knot[i];
  for(k=0; k<top; k++)
    knot[k] = knot[k] - newOrigin;
  
  /* M_PI and M_PI_2 are constants in math.h */
  /* Rotate so that vertex j is on the x-axis */
  double theta = -atan( knot[j][1] / knot[j][0] );
  if( knot[j][0] < 0 )
    theta = -M_PI + theta;

  double phi = (knot[j][2] / ( sqrt( knot[j][0]*knot[j][0] + 
	       		  knot[j][1]*knot[j][1] + knot[j][2]*knot[j][2] )) );
  if(phi >= 1)
    phi = 0;
  else
    phi = acos(phi);

  /* phi correction gives angle to x-axis */
  phi = M_PI_2 - phi;
  knot.rotatexy(theta);
  knot.rotatezx(phi);
}  

/* translates knot so that v1 is at the origin and direction lies
   on the positive x-axis (y and z coords are 0)*/
void PointsNormalize(knots & knot, const point & v1, const point & direction)
{
  int k, top = knot.vnum();

  /* Translate */
  point newOrigin = v1;
  for(k=0; k<top; k++)
    knot[k] = knot[k] - newOrigin;
  
  point v2 = direction - newOrigin;

  /* M_PI and M_PI_2 are constants in math.h */
  /* Rotate so that vertex j is on the x-axis */
  double theta = -atan( v2[1] / v2[0] );
  if( v2[0] < 0 )
    theta = -M_PI + theta;

  double phi = (v2[2] / ( sqrt( v2[0]*v2[0] +  v2[1]*v2[1] + v2[2]*v2[2] )) );
  if(phi >= 1)
    phi = 0;
  else
    phi = acos(phi);

  /* phi correction gives angle to x-axis */
  phi = M_PI_2 - phi;
  knot.rotatexy(theta);
  knot.rotatezx(phi);
}  

/* does the rotation of the part of the knot between vertex i and j */
void Crank(knots & knot, const double & angle, const int & i, const int & j)
{
  knot.rotateyz(angle,i,j);
}

/* computes maximum angle that guarantees that the knot type is preserved */
double ComputeMaxAngle(const knots & knot, const int & i, const int & j)
{
  double d = knot[j][0];
  if (d<=0)
    {
      cerr << "whoops, d = 0!\n";
      Save(knot, "ERROR");
      exit(0);
    }
      
  double b;
  if (2*(j-i) > knot.vnum())
    b = knot.lenvtr[0] * double(knot.vnum() - j + i);
  else
    b = knot.lenvtr[0] * double(j-i);

  b = b*b;
  d = d*d;
  double theta = ((2*knot.injrad()) / (sqrt(b + d)));
  /* make sure that we aren't taking asin of a big number */
  if(theta >= 1)
    {
      theta = M_PI_2;
    }
  else
    theta = 2 * asin(theta);

  return theta;
}

/* equilateral version of run */
void EqRun(knots & knot, int runs)
{
  double chance;
  int i, j;

  for(; runs > 0; runs--)
    {
      knots tempknot(knot);
      double tempRopeLength = knot.RopeLength();
      /* choose random non-consecutive vertices */
      do 
	{
	  i = int(drand48()*knot.vnum());
	  j = int(drand48()*knot.vnum());

	  assert (j < knot.vnum());

	  /* reorder if necessary */
	  if(j < i)  
	    {
	      int k = i;
	      i = j;
	      j = k;
	    }
	  /* make sure that we have a legal pairing of i and j */
	} while( (i == j) || ((j - i) == 1) || 
		 (i == 0 && j == (knot.vnum()-1))); 

      Normalize(knot,i,j);  /* translate and rotate */

      /* the 2 is to quicken the following computation */
      double maxangle = 2 * ComputeMaxAngle(knot, i, j);

      chance = ((drand48() - .5) * maxangle);

      Crank(knot,chance,i,j);
      knot.Compute();

      /* compareThickness(knot, tempknot); */
      if (knot.RopeLength() > tempRopeLength)
	knot = tempknot;
    }

  cout << "Final knot data:" << endl;
  knot.Display();
}

/* equilateral version of run */
void Tidy(knots & knot, int runs)
{
  double chance;
  int i, top = knot.vnum();
  
  for(; runs > 0; runs--)
    {
      knots tempknot(knot);
      double tempRopeLength = knot.RopeLength();

      for(i=0; i<top; i++)
	{
	  Normalize(knot,i,i+2);  /* translate and rotate */
	  
	  chance = ((drand48() - .5) * .314);
	  
	  Crank(knot,chance,i,i+2);
	  knot.Compute();
	  
	  /* compareThickness(knot, tempknot); */
	  if (knot.RopeLength() > tempRopeLength)
	    knot = tempknot;
	}
    }
  cout << "Final knot data:" << endl;
  knot.Display();
}

/* equilateral version of conjugate */
void EqConjugateRun(knots & knot, int runs)
{
  double chance;
  int i, j;

  for(; runs > 0; runs--)
    {
      double tempRopeLength;
      knots tempknot;
      /* choose random non-consecutive vertices */
      do 
	{
	  i = int(drand48()*knot.vnum());
	  j = int(drand48()*knot.vnum());

	  /* reorder if necessary */
	  if(j < i)  
	    {
	      int k = i;
	      i = j;
	      j = k;
	    }

	} while( (i == j) || ((j - i) == 1) || 
		 (i == 0 && j == (knot.vnum()-1))); 

      Normalize(knot,i,j);  /* translate and rotate */
      double maxangle = 2 * ComputeMaxAngle(knot, i, j);

      /* the 2 is to quicken the following computation */
      chance = ((drand48() - .5) * maxangle);

      /* compare the ropelengths and keep trying the same vertices as long
	 as it keeps working */
      do{
	tempRopeLength = knot.RopeLength();
	tempknot = knot;
	Crank(knot,chance, i,j);
	knot.Compute();
      } while (knot.RopeLength() <= tempRopeLength);

      knot = tempknot;
    }

  cout << "Final knot data:" << endl;
  knot.Display();
}

/* equilateral version of persistant, can be slow */
void EqPersistentRun(knots & knot, int runs)
{
  double chance;
  int i, j;

  for(; runs > 0; runs--)
    {
      double tempRopeLength = knot.RopeLength();
      knots tempknot(knot);
      /* choose random non-consecutive vertices */
      do 
	{
	  i = int(drand48()*knot.vnum());
	  j = int(drand48()*knot.vnum());

	  /* reorder if necessary */
	  if(j < i)  
	    {
	      int k = i;
	      i = j;
	      j = k;
	    }

	} while( (i == j) || ((j - i) == 1) || 
		 (i == 0 && j == (knot.vnum()-1))); 

      Normalize(knot,i,j);  /* translate and rotate */
      double maxangle = 2 * ComputeMaxAngle(knot, i, j);

      /* the 2 is to quicken the following computation */
      chance = ((drand48() - .5) * maxangle);

      Crank(knot,chance, i,j);
      knot.Compute();

      /* compare the ropelengths and keep trying the same vertices as long
	 as it keeps working, if it doesn't work then try fractions */
      if(knot.RopeLength() <= tempRopeLength)
	{ 
	  do{
	    tempRopeLength = knot.RopeLength();
	    tempknot = knot;
	    Crank(knot, chance, i, j);
	    knot.Compute();
	  } while (knot.RopeLength() <= tempRopeLength);
	}
      else
	{
	  int h=0;
	  do{
	    h++;
	    knot = tempknot;
	    chance *= .5;
	    Crank(knot, chance, i, j);
	    knot.Compute();
	  } while( (knot.RopeLength() > tempRopeLength) && (h<5) );
	}	  

      knot = tempknot;
    }

  cout << "Final knot data:" << endl;
  knot.Display();
}

/* equilateral version of annealing */
void EqAnnealingRun(knots & knot, double temp, 
		  const double & dissipation, int runs)
{
  double chance;
  int i, j;
  knots shortest(knot);
  double shortestLength = knot.RopeLength();
  double qtemp = 3000/temp;   /* for quickness */
  int attempts=0,anneals=0;

  for(; runs > 0; runs--)
    {
      knots tempknot(knot);
      double tempRopeLength = knot.RopeLength();
      /* choose random non-consecutive vertices */
      do 
	{
	  i = int(drand48()*knot.vnum());
	  j = int(drand48()*knot.vnum());

	  assert (j < knot.vnum());

	  /* reorder if necessary */
	  if(j < i)  
	    {
	      int k = i;
	      i = j;
	      j = k;
	    }
	} while( (i == j) || ((j - i) == 1) || 
		 (i == 0 && j == (knot.vnum()-1))); 

      Normalize(knot,i,j);  /* translate and rotate */
      double maxangle = 2 * ComputeMaxAngle(knot, i, j);

      /* the 2 is to quicken the following computation */
      chance = ((drand48() - .5) * maxangle);

      Crank(knot,chance, i,j);
      knot.Compute();
      /* compareThickness(knot, tempknot); */
      if (knot.RopeLength() > tempRopeLength)
	{
	  attempts++;
	  double probability = drand48();
	  /* percent is percent change from old rope length and is negative */
	  double percent = (tempRopeLength-knot.RopeLength())/tempRopeLength;
	  if(exp(percent*qtemp)<drand48())
	    {
	      anneals++;
	      knot = tempknot;
	    }
	}
      else if(knot.RopeLength() <= shortestLength)
	{
	  shortestLength = knot.RopeLength();
	  shortest = knot;
	}
    }

  string filename;

  cout << "Save shortest to file: ";
  cin >> filename;

  Save(shortest, filename);
  cout << "Final knot data:" << endl;
  cout << "Shortest was " << shortest.RopeLength() << endl;
  cout << "Attempts = " << attempts << "  Anneals = " << anneals
       << "  % = " << 1-double(anneals)/attempts << endl;
  knot.Display();
}


/* --------------------------------------------------------------------- */
/* Input and Output of Files                                             */
/* --------------------------------------------------------------------- */

void Load(knots & knot, string & file)
{
  ifstream fin(file.c_str());
  while(! (fin.good()))
    {
      cout << "File name is invalid!\nEnter new file name: ";
      cin >> file;
      fin.close();
      fin.open(file.c_str());
    }

  knots temp(knot);

  if(fin >> knot)
    { 
      cout << file << " with " << knot.vnum() << " vertices loaded.\n";
    }
  else
    {
      knot = temp;
      cout << "Bad knot file!  Try another file." << endl;
    }

}

void LoadMillett(knots & knot, const string & file)
{
  cout << "Not functional!\n";
}

void Save(const knots & knot, const string & file)
{
  ofstream fout(file.c_str());

  fout << knot;

  cout << "Knot successfully saved to " << file << endl;
}

void SaveMillett(const knots & knot, const string & file)
{
  int top = knot.vnum()-1;
  
  ofstream fout(file.c_str());
  fout.setf(ios::fixed);
  fout.precision(10);

  fout << "{";
  for(int i=0; i<top; i++)
    {
      fout << "{" << knot[i][0] << ", " << knot[i][1]
	<< ", " << knot[i][2] << "},\n";
    }

      fout << "{" << knot[top][0] << ", " << knot[top][1]
	<< ", " << knot[top][2] << "}}\n";  
}

/* doubles the number of sides by bisecting each side */
void DoubleSides(knots & knot)
{
  knots old(knot);

  int top = knot.vnum();

  knot.setSize(top*2);

  for(int i=0; i<top; i++)
    {
      knot[2*i] = old[i];
      knot[2*i+1] = .5*old[i]+.5*old[((i+1)%top)];
    }

  knot.Compute();
}


/* This doesn't work as is */
knots MakeEquilateral(const knots & knot)
{
  double targetLength=0;

  while(targetLength <= 0)
    {
      cout << "Enter the target side length: ";
      cin >> targetLength;
    }

  knots temp(knot);

  int top = (knot.vnum() - 1);

  for(int i=0; i<top; i++)
    {
      Normalize2(temp,i,i+1);
      temp[i+1] = point(targetLength,0,0);
    }

  Normalize2(temp,0,(top-1));

  temp.Compute();
  cout << temp << endl;
  temp.SideData();

  return temp;
}

knots MakeEquilateral2(knots & knot)
{
  double targetLength=0;
  int top = (knot.vnum() - 1);
  int nexttolast = top - 1;

  cout.precision(16);

  while(targetLength <= 0)
    {
      cout << "Enter the target edge length: ";
      cin >> targetLength;
    }

  cout << "original: " << knot << endl;

  Normalize2(knot,0,1);
  knot.Compute();

  cout << "normalized: " << knot << endl;

  /* 0th should be at origin and 1st should be on x-axis */
  /* Now I want the nexttolast on first quadrant of xy plane */
  cout << "0th at origin, 1st on x-axis\n";
  cout << "0th: " << knot[0] << "\n1st: " << knot[1] << endl;

  double b = knot[nexttolast][1], c = knot[nexttolast][2];
  
  double theta = b/(sqrt(b*b+c*c));
  if(theta <= -1)
    theta = -1;
  else if(theta >= 1)
    theta = 1;
  else
    theta = theta;

  if(c < 0)
    knot.rotaterawnegyz(theta);
  else
    knot.rotaterawyz(theta);

  cout << "Also, nexttolast should have z=0\n";
  cout << "0th: " << knot[0] << "\n1st: " << knot[1] << "\nnexttolast: " 
       << knot[nexttolast] << endl;

  knot.Compute();

  cout << "rotated: " << knot << endl;

  knots temp(knot);

  temp[0] = point(0.0,0.0,0.0);
  for(int i=1; i<top; i++){
    temp[i] = temp[i-1] + (targetLength/knot.lenvtr[i-1])*knot.side[i-1];
  }

  /* Now rotate outnormed to the right spot */
  /* First renormalize temp so that nexttolast is on positive xy plane */
  b = temp[nexttolast][1], c = temp[nexttolast][2];

  theta = b/(sqrt(b*b+c*c));
  if(theta <= -1)
    theta = 0;
  else if(theta >= 1)
    theta = -M_PI;
  else
    theta = acos(theta);

  if(c < 0)
    theta = -theta;

  temp.rotateyz(theta);
  cout << "Also, nexttolast should have z=0\n";
  cout << "0th: " << temp[0] << "\n1st: " << temp[1] << "\nnexttolast: " 
       << temp[nexttolast] << endl;

  double gap = (temp[nexttolast]-temp[0]).norm();

  cout << "Last two are: " << gap << " apart.\n";
  cout << "This should be > 0 and < " << 2*targetLength << endl;

  if(gap == 0){
    cerr << "The gap is exactly 0, that's bad!\n";
    return knot;
  }
  else if(gap >= 2*targetLength){
    cerr << "The gap is > " << 2*targetLength << ", and that's terrible!\n";
    return knot;
  }

  double needed = sqrt(targetLength*targetLength - 0.25*gap*gap);
  cout << "needed length: " << needed << endl;

  point knotmidpoint = 0.5*(knot[0] + knot[nexttolast]);
  point out = knot[top] - knotmidpoint;
  point outnormed = (needed/out.norm())*out;

  cout << "length of outnormed: " << outnormed.norm() << endl;

  point tempmidpoint = 0.5*(temp[0] + temp[nexttolast]);
  temp[top] = tempmidpoint + outnormed;

  temp.Compute();
  cout << "New:\n";
  cout << temp << endl;
  temp.SideData();
  cout << "Old:\n";
  cout << knot << endl;
  knot.SideData();

  return temp;
}

      
void Normalize2(knots & knot, const int & i, const int & j)
{
  int k, top = knot.vnum();

  /* Translate */
  point newOrigin = knot[i];
  for(k=0; k<top; k++)
    knot[k] = knot[k] - newOrigin;
  // long double a = knot[j][0], b = knot[j][1], c = knot[j][2];
  double a = knot[j][0], b = knot[j][1], c = knot[j][2];


  /* M_PI and M_PI_2 are constants in math.h */
  /* Rotate so that vertex j is on the x-axis */
  // long double theta = -atan(b/a);
  double theta = -atan(b/a);
  if( a < 0 )
    theta = -M_PI + theta;

  // long double phi = (c / ( sqrt(a*a + b*b + c*c)));
  double phi = (c / ( sqrt(a*a + b*b + c*c)));
  if(phi >= 1)
    phi = 0;
  else
    phi = acos(phi);


  /* phi correction gives angle to x-axis */
  phi = M_PI_2 - phi;

  // cout << "theta " << theta << " phi " << phi << endl;

  knot.rotatexy(theta);
  knot.rotatezx(phi);
  knot[j][1] = knot[j][2] = 0;
}

/* show information on vertex including point and angles and sides adjacent */
void ShowVertexInfo(const knots & knot)
{
  cout << "Examine vertex number: ";
  int number;
  cin >> number;

  if( (number > 0) && (number < knot.vnum()) )
    {
      cout.setf(ios::fixed);
      cout.precision(14);

      cout << "Vertex #" << number << endl;
      cout << "vertex[" << number << "] = " << knot[number] << endl;
      cout << "side[" << (number-1+knot.vnum()) % knot.vnum() 
	<< "] = " << knot.side[number-1] << endl;
      cout << "side[" << number
	<< "] = " << knot.side[number] << endl;
      cout << "difference is " << knot.side[number]-knot.side[number-1]
	<< endl;
      cout << "angle[" << number << "] = " << knot.angle[number] << endl;
    }
  else
    cout << "Bad vertex number!\n";
}


/* resize knot */
void Resize(knots & knot, const double & newsize)
{
  double ratio = (newsize / knot.Length());
  int top = knot.vnum();

  for(int i=0; i<top; i++)
    knot[i] = ratio*knot[i];
  
  knot.Compute();
}


void Standardize2(knots & knot, double & d1, double & d2, double & d3)
{
  /* this function finds the bounds and volume of an ellipsoid 
     which contains the knot */
  double temp;
  int i,j,l1,l2,w,h, top = knot.vnum();
  d1 = d2 = d3 = 0;

  Resize(knot,100);

  for(i=0; i<top; i++)
    for(j=(i+1); j< top; j++)
      {
	temp = (knot[i]-knot[j]).norm();
	if(temp > d1)
	  {
	    d1 = temp;
	    l1 = i;
	    l2 = j;
	  }
      }

  Normalize(knot, l1, l2);

  temp = knot[l2].norm();
  temp = -.5*temp;

  for(i=0; i<top; i++)
    knot[i][0] = knot[i][0] + temp;

  /* determine vertex that is second furthest from origin now */
  for(i=0; i<top; i++)
    {
	temp = (knot[i]).norm();
	if( (temp > d2) && (i!=l1) && (i!=l2) )
	  {
	    d2 = temp;
	    w = i;
	  }
      }  

  /* project onto yz-plane and find angle with y-axis */
  point p = knot[w];
  p[0] = 0;
  double theta = acos(p[1]/(p.norm()));
  if(p[2]>0)
    theta = -theta;

  knot.rotateyz(theta,-1,top);

  knots projection(knot);

  /* projected onto z-axis */
  for(i=0; i<top; i++)
    projection[i][0] = projection[i][1] = 0;

  for(i=0; i<top; i++)
    {
      temp = fabs(projection[i][2]);
      if(temp > d3)
	{
	  d3 = temp;
	  h = i;
	}
    }  

  knot.Compute();
  
  /* determine the eccentricity of an ellipsoid engulfing K 
     (a,b) is (x,y) coords of knot[w] after projected to xy-plane             
     then (a,b,c) is (x,y,z) coords of knot[h]
     looking for x^2/r^2+y^2/s^2+z^2/t^2=1 */

  double a,b,c,r,s,t;
  a = knot[w][0];
  b = knot[w][1];
  cout << "(a,b) = (" << a << " , " << b << ")\n"; 
  r = knot[l2][0];
  cout << "r = " << r << endl;
  s = sqrt((b*b*r*r)/(r*r-a*a));

  a = knot[h][0];
  b = knot[h][1];
  c = knot[h][2];
  cout << "(a,b,c) = (" << a << " , " << b << " , " << c << ")\n"; 

  t = sqrt((c*c*r*r*s*s)/(r*r*s*s - s*s*a*a - r*r*b*b));

  cout << "r = " << r << endl;
  cout << "s = " << s << endl;
  cout << "t = " << t << endl;
  cout << "Volume = " << 4*r*s*t*M_PI/3 << endl;

  cout << "Eccentricity " << d1 << " " << d2 << " " << d3 << endl;
  cout << "sides: " << l1 << " " << l2 << endl;
  cout << "sides: " << w << endl;
  cout << "sides: " << h << endl;
}

void Standardize(knots & knot, double & d1, double & d2, double & d3)
{
  /* Finds the box dimensions and the volume and surface area of the box */

  double temp;
  int i,j,v1,v2,v3,v4,v5,v6, top = knot.vnum();
  d1 = d2 = d3 = 0;

  Resize(knot,100);

  for(i=0; i<top; i++)
    for(j=(i+1); j< top; j++)
      {
	temp = (knot[i]-knot[j]).norm();
	if(temp > d1)
	  {
	    d1 = temp;
	    v1 = i;
	    v2 = j;
	  }
      }

  Normalize(knot, v1, v2);

  /* Save(knot,"lengthknot");
     cout << "v1 " << v1 << " v2 " << v2 << endl; */

  /* project onto yz-plane and */
  /* determine vertex pair that is second furthest apart */
  knots projection(knot);
  for(i=0; i<top; i++)
    projection[i][0] = 0;

  for(i=0; i<top; i++)
    for(j=i+1; j<top; j++)
      {
	temp = (projection[i]-projection[j]).norm();
	if( (temp > d2) && ((i!=v1) || (i!=v2)) )
	  {
	    d2 = temp;
	    v3 = i;
	    v4 = j;
	  }
      }  

  /* Save(projection,"widthknot");
     cout << "v3 " << v3 << " v4 " << v4 << endl; */

  /* rotate so that p is parallel to the xy-plane*/
  point p = knot[v4]-knot[v3];
  p[0] = 0;
  double theta = acos(p[1]/(p.norm()));
  if(p[2]>0)
    {
      /* cout << theta << " " << p << endl; */
      theta = -theta;
      /* cout << "it's negative\n"; */
    }
  knot.rotateyz(theta,-1,top);
  knot.Compute();
  /* cout << "side 1 is " << knot[v2]-knot[v1] << endl;
     cout << "side 2 is " << knot[v4]-knot[v3] << endl; */

  /* projected onto z-axis */
  projection = knot;
  for(i=0; i<top; i++)
    projection[i][1] = projection[i][0] = 0;

  for(i=0; i<top; i++)
    for(j=i+1; j<top; j++)
      {
	temp = (projection[i]-projection[j]).norm();
	if( (temp > d3) && ((i!=v1) || (j!=v2)) && ((i!=v3) || (j!=v4)) )
	  {
	    d3 = temp;
	    v5 = i;
	    v6 = j;
	  }
      }  

  /* Save(projection,"heightknot");
     cout << "v5 " << v5 << " v6 " << v6 << endl; */

  knot.Compute();
  
  cout << "Box dimensions: " << d1 << " X " << d2 << " X " << d3 << endl;
  cout << "Box volume: " << d1*d2*d3 << endl;
  cout << "Box surface area: " << 2*d1*d2 + 2*d2*d3 + 2*d1*d3 << endl;
  cout << "sides: " << v1 << " " << v2 << endl;
  cout << "sides: " << v3 << " " << v4 << endl;
  cout << "sides: " << v5 << " " << v6 << endl;
}


void ACWN(const knots & knot, double & acn, double & awn)
{
  acn = awn = 0;
  point q[4];
  int top = knot.vnum();
  double c[4];
  double sum;
  double denominator;
  int ip1, jp1;

  for(int i=0; i<top; i++)
    for(int j=i+2; j<top; j++)
      {
	if( (i!=0) || (j!=(top-1)) )
	  {
	    ip1 = ((i+1)%top);
	    jp1 = ((j+1)%top);

	    // cout << i << " " << j << endl;
	    // cout << ip1 << " " << jp1 << endl;

	    q[0] = knot[j] - knot[ip1];
	    q[1] = knot[jp1] - knot[ip1];
	    q[2] = knot[jp1] - knot[i];
	    q[3] = knot[j] - knot[i];
	    
	    denominator = ((q[3] & q[0]).norm() * (q[0] & q[1]).norm());
	    if(denominator != 0){
	      c[0] = ( ((q[3] & q[0]) * (q[0] & q[1])) / (denominator));
	    }
	    else
	      c[0] = 0;
	    // cout << "301: " << denominator << endl;
	    // cout << ((q[3] & q[0]) * (q[0] & q[1])) << endl;
	    // cout << c[0] << endl;

	    denominator = ((q[0] & q[1]).norm() * (q[1] & q[2]).norm());
	    if(denominator != 0){
	      c[1] = ( ((q[0] & q[1]) * (q[1] & q[2])) / (denominator));
	    }
	    else
	      c[1] = 0;
	    // cout << "012: " << denominator << endl;
	    // cout << ((q[0] & q[1]) * (q[1] & q[2])) << endl;
	    // cout << c[1] << endl;
	    
	    denominator = ((q[1] & q[2]).norm() * (q[2] & q[3]).norm());
	    if(denominator != 0){
	      c[2] = ( ((q[1] & q[2]) * (q[2] & q[3])) / (denominator));
	    }
	    else
	      c[2] = 0;
	    // cout << "123: " << denominator << endl;
	    // cout << ((q[1] & q[2]) * (q[2] & q[3])) << endl;
	    // cout << c[2] << endl;
	    
	    denominator = ((q[2] & q[3]).norm() * (q[3] & q[0]).norm());
	    if(denominator != 0){
	      c[3] = ( ((q[2] & q[3]) * (q[3] & q[0])) / (denominator));
	    }
	    else
	      c[3] = 0;
	    // cout << "230: " << denominator << endl;
	    // cout << ((q[2] & q[3]) * (q[3] & q[0])) << endl;
	    // cout << c[3] << endl;
	    

	    for(int k=0; k<4; k++)
	      {
		if (c[k] > 1)
		  c[k] = 0;
		else if (c[k] < -1)
		  c[k] = M_PI;
		else
		  c[k] = acos(c[k]);
	      }
	    
	    sum = c[0] + c[1] + c[2] + c[3];
	    sum = 1 - sum/(2*M_PI);

	    if(sum < 0.0)
	      sum = 0.0;
	    else if(sum > 1.0)
	      sum = 1.0;

	    // cout << "i " << i << " j " << j << " sum " << sum << endl << endl;

	    acn += sum;
	    if( ((q[0] & q[1]) * q[2]) < 0.0)
	      awn -= sum;
	    else
	      awn += sum;
	  }
      }
}


// void ACWN(const knots & knot, double & acn, double & awn)
// {
//   acn = awn = 0;
//   point q[4];
//   int top = knot.vnum();
//   double c[4];
//   double sum;
//   double denominator;    

//   for(int i=0; i<top; i++)
//     for(int j=i+2; j<top; j++)
//       {
// 	cout << "in loop: i " << i << "  j " << j << endl;

// 	if( (i!=0) || (j!=(top-1)) )
// 	  {
// 	    cout << "okay: i " << i << "  j " << j << endl;

// 	    q[0] = knot[(j+1)%top] - knot[i];
// 	    q[1] = knot[j%top] - knot[i];
// 	    q[2] = knot[j%top] - knot[(i+1)%top];
// 	    q[3] = knot[(j+1)%top] - knot[(i+1)%top];
	    
// 	    denominator = ((q[3] & q[0]).norm() * (q[0] & q[1]).norm());
// 	    if(denominator > .000000000000001)
// 	      c[0] = ( ((q[3] & q[0]) * (q[0] & q[1])) / (denominator));
// 	    else
// 	      c[0] = 0;
	     
// 	    denominator = ((q[0] & q[1]).norm() * (q[1] & q[2]).norm());
// 	    if(denominator > .000000000000001)
// 	      c[1] = ( ((q[0] & q[1]) * (q[1] & q[2])) / (denominator));
// 	    else
// 	      c[1] = 0;

// 	    denominator = ((q[1] & q[2]).norm() * (q[2] & q[3]).norm());
// 	    if(denominator > .000000000000001)
// 	      c[2] = ( ((q[1] & q[2]) * (q[2] & q[3])) / (denominator));
// 	    else
// 	      c[2] = 0;

// 	    denominator = ((q[2] & q[3]).norm() * (q[3] & q[0]).norm());
// 	    if(denominator > .000000000000001)
// 	      c[3] = ( ((q[2] & q[3]) * (q[3] & q[0])) / (denominator));
// 	    else
// 	      c[3] = 0;

// 	    for(int k=0; k<4; k++)
// 	      {
// 		if (c[k] > 1)
// 		  c[k] = 1;
// 		if (c[k] < -1)
// 		  c[k] = -1;
// 	      }
	    
// 	    sum = acos(c[0]) + acos(c[1]) + acos(c[2]) + acos(c[3]);
// 	    sum = 1 - sum/(2*M_PI);
	    
// 	    acn += sum;
// 	    if( ((q[0] & q[1]) * q[2]) < 0)
// 	      awn -= sum;
// 	    else
// 	      awn += sum;
// 	  }
//       }
// }


void AcnAnneal(knots & knot, const double & temperature, int runs)
{
  double chance;
  int i, j;
  knots leastACN(knot);
  double tempACN, dummy,ACN;
  ACWN(knot,tempACN,dummy);
  ACN = tempACN;
  knots tempknot(knot);
  double lowestACN = ACN;
  double number = 3000/temperature;   /* to quicken */

  for(; runs > 0; runs--)
    {
      /* choose random non-consecutive vertices */
      do 
	{
	  i = int(drand48()*knot.vnum());
	  j = int(drand48()*knot.vnum());

	  /* reorder if necessary */
	  if(j < i)  
	    {
	      int k = i;
	      i = j;
	      j = k;
	    }
	} while( (i == j) || ((j - i) == 1) || 
		 (i == 0 && j == (knot.vnum()-1))); 

      Normalize(knot,i,j);  /* translate and rotate */
      double maxangle = 2 * ComputeMaxAngle(knot, i, j);

      /* the 2 is to quicken the following computation */
      chance = ((drand48() - .5) * maxangle);

      Crank(knot,chance, i,j);
      knot.Compute();
      ACWN(knot,ACN,dummy);
      /* compareThickness(knot, tempknot); */
      if (ACN > tempACN)
	{
	  double probability = drand48();
	  /* percent = percent change from old rope length and is negative */
	  double percent = (tempACN-ACN)/tempACN;
	  if(exp(percent*number)<drand48())
	    {
	      knot = tempknot; /* back to the original */
	    }
	  else
	    {
	      tempACN = ACN;
	      tempknot = knot;
	    }
	}
      else if(ACN < lowestACN)
	{
	  lowestACN = tempACN = ACN;
	  leastACN = knot;
	  tempknot = knot;
	}
      else
	{
	  tempACN = ACN;
	  tempknot = knot;
	}
    }
  
  string filename;
  
  cout << "Save knot with lowest ACN to file: ";
  cin >> filename;
  
  Save(leastACN, filename);
  cout << "Final knot data:" << endl;
  cout << "Lowest ACN was " << lowestACN << endl;
  knot.Display();
}

void HowThickSmooth(const knots & knot)
{
  double lengthSaved = 0;  /* at each vertex, save R*(2*tan(angle/2)-angle) */
  double biggestAngle = 0;  
  double smoothLength = 0;  /* knot.Length() - lengthSaved */
  double injrad = knot.injrad();
  double angle;  /* temporary storage for the angle */
  int top = knot.vnum();

  for(int k=0; k<top; k++){
    angle = knot.angle[k];
    lengthSaved += 2*tan(0.5*angle) - angle;
    if(angle > biggestAngle)
      biggestAngle = angle;
  }

  lengthSaved *= injrad;

  double smoothinjrad = injrad*(2 - 1/(cos(biggestAngle/2)));
  smoothLength = knot.Length() - lengthSaved;

  cout << "Maximum angle = " << biggestAngle << endl;
  cout << "Polygon Length = " << knot.Length() << endl;
  cout << "Length saved on smooth = " << lengthSaved << endl;
  cout << "Smooth knot's length = " << smoothLength << endl;
  cout << "Poly injectivity radius = " << injrad << endl;
  cout << "Smooth injectivity radius = " << smoothinjrad << endl;

  cout << "Ropelength upper bound = " << smoothLength / smoothinjrad << endl;
}

void HowThickSmooth2(const knots & knot)
{
  double lengthSaved = 0;
  double biggestAngle = 0;
  double halfbiggestAngle;
  int top = knot.vnum();
  double injrad = knot.injrad();
  double length = knot.Length();
  double newlength;
  double minrad = knot.MinRad();

  for(int k=0; k<top; k++)
    {
      if(knot.angle[k] > biggestAngle)
	biggestAngle = knot.angle[k];
      lengthSaved += ((2*tan(knot.angle[k]*.5)) - knot.angle[k]);
    }
  lengthSaved *= injrad;
  newlength = length - lengthSaved;

  halfbiggestAngle = .5*biggestAngle;

  cout << "MinRad = " << minrad << endl;
  cout << "Maximum angle = " << biggestAngle << endl << endl;
  cout << "Polygon Length = " << length << endl;
  cout << "Length saved on smooth = " << lengthSaved << endl;
  cout << "Smooth knot's length = " << newlength << endl << endl;

  double cornerError = injrad*((1/cos(halfbiggestAngle)) - 1);
  double withcornerir = injrad - cornerError;

  cout << "Injectivity radius = " << injrad << endl;
  cout << "Corner error = " << cornerError << endl;
  cout << "dc with corner error = " << withcornerir << endl;

  if(withcornerir <= 0)
    cout << "Approximation Theorem tells us nothing.\n";
  else
    {
      cout << newlength/injrad << " <= Smooth ropelength <= "
	   << newlength/withcornerir << endl;
    }
}

void HowThickSmoothFactor(const knots & knot)
{
  cout << "Factor (0<f<1): ";
  double f;
  cin >> f;
  // knot.ChangeMinRadFactor(f);

  double invf = 1/f;

  double lengthSaved = 0;
  double biggestAngle = 0;
  double halfbiggestAngle;
  int top = knot.vnum();
  double injrad = knot.injrad();
  double length = knot.Length();
  double newlength;
  double minrad = knot.MinRad();

  for(int k=0; k<top; k++)
    {
      if(knot.angle[k] > biggestAngle)
	biggestAngle = knot.angle[k];
      lengthSaved += ((2*tan(knot.angle[k]*.5)) - knot.angle[k]);
    }
  lengthSaved *= injrad;
  newlength = length - lengthSaved;

  halfbiggestAngle = .5*biggestAngle;

  cout << "MinRad = " << minrad << endl;
  cout << "Maximum angle = " << biggestAngle << endl << endl;
  cout << "Polygon Length = " << length << endl;
  cout << "Length saved on smooth = " << lengthSaved << endl;
  cout << "Smooth knot's length = " << newlength << endl << endl;

  double cornerError = 2*injrad*((1/cos(halfbiggestAngle)) - 1);
  double withcornerir = injrad - cornerError;

  cout << "Injectivity radius = " << injrad << endl;
  cout << "Corner error = " << cornerError << endl;
  cout << "dc with corner error = " << withcornerir << endl;

  if(withcornerir <= 0)
    cout << "Approximation Theorem tells us nothing.\n";
  else
    {
      cout << newlength/injrad << " <= Smooth ropelength <= "
	   << newlength/withcornerir << endl;
    }
}






void CheckCriticality(const knots & knot)
{
  int top = knot.vnum(), i, j;
  double d;
  double smallestCrit = 65000, biggestCrit = 0;
  vectors cdist(top,65000);
  double sumangle=0;
  int start, end;

  for(i=0; i<top; i++)
    {
      sumangle = knot.angle[(i+1)%top];
      start = i+1;
      while(sumangle < M_PI)
	{
	  start++;
	  if(start >= top)
	    start %= top;
	  sumangle += knot.angle[start];
	}
      start ++;
      
      sumangle = knot.angle[i];
      end = i;
      while(sumangle < M_PI)
	{
	  end --;
	  if(end < 0)
	    end += top;
	  sumangle += knot.angle[end];
	}

      if(end < 0)
	end += top;

      if(start < end)
	{
	  for(j=start; j<end; j++)
	    {
	      if( (i!=j) && (((i+1)%top) != j) && (i !=((j+1)%top))){
		NewFindRate(knot,i,j,d);
		if(d<cdist[i])
		  cdist[i] = d;
	      }
	    }
	}
      else if(end < start)
	{
	  for(j=start; j<top; j++)
	    {
	      if( (i!=j) && (((i+1)%top) != j) && (i !=((j+1)%top))){
		NewFindRate(knot,i,j,d);
		if(d<cdist[i])
		  cdist[i] = d;
	      }
	    }  

	  for(j=0; j<end; j++)
	    {
	      if( (i!=j) && (((i+1)%top) != j) && (i !=((j+1)%top))){
		NewFindRate(knot,i,j,d);
		if(d<cdist[i])
		  cdist[i] = d;
	      }
	    }
	}
      else
	{
	  cerr << "May be problems with total curvature stuff!\n";
	  d = (knot[i] - knot[start]).norm();
	  if(d<cdist[i])
	    cdist[i] = d;
	}
    }

  smallestCrit = 65000;
  biggestCrit = 0;

  for(i=0; i<top; i++)
    {
      if(cdist[i] < smallestCrit)
	smallestCrit = cdist[i];
      if(cdist[i] > biggestCrit)
	biggestCrit = cdist[i];
    }
  
  printf("%s\t%s\t%s\t\t%s\n", "Vertex", "    Dist", "    r_-", "    r_+");

  for(i=0; i<top; i++)
    {
      printf("%i\t%14.8f\t%14.8f\t%14.8f\n", i, cdist[i], 
	   knot.lenvtr[((i-1+top)%top)]/(2*tan(knot.angle[i]*.5)),
	   knot.lenvtr[i]/(2*tan(knot.angle[i]*.5)));
    }

  cout << "Smallest = " << smallestCrit << endl;
  cout << "Biggest = " << biggestCrit << endl;
}



void CheckCriticality2(const knots & knot)
{
  int top = knot.vnum(), i, j;
  double d;
  double smallestCrit = 65000, biggestCrit = 0;
  vectors cdist(top,65000);

  for(i=0; i<top; i++)
    {
      for(j=0; j<top; j++)
	{
	  if( (i!=j) && (((i+1)%top)!=j) && (((j+1)%top) != i) )
	    {
	      if(NewFindRate(knot,i,j,d) && (d<cdist[i]))
		cdist[i] = d;
	    
	      point e= knot[i]-knot[j];
	      double a=(e*knot.side[j])/(knot.side[j]*knot.side[j]);
	      if (0<=a && a<=1)
		{
		  if (IsDC(knot,i,j,a))
		    {
		      double d=(e-(a*knot.side[j])).norm();
		      if ( d < cdist[i] )
			cdist[i] = d;
		    }
		}
		
	      e= knot[j]-knot[i];
	      a=(e*knot.side[i])/(knot.side[i]*knot.side[i]);
	      if (0<=a && a<=1)
		{
		  if (IsDC(knot,j,i,a))
		    {
		      double d=(e-(a*knot.side[i])).norm();
		      if ( d < cdist[i] )
			cdist[i] = d;
		    }
		}

	      if (IsDC(knot,i,j,0) && IsDC(knot,j,i,0))
		{
		  d=(knot[i]-knot[j]).norm();
		  if ( d < cdist[i] )
		    cdist[i] = d;
		}
	    }
	}

    }

  smallestCrit = 65000;
  biggestCrit = 0;

  for(i=0; i<top; i++)
    {
      if(cdist[i] < smallestCrit)
	smallestCrit = cdist[i];
      if(cdist[i] > biggestCrit)
	biggestCrit = cdist[i];
    }
  
  printf("%s\t%s\t%s\t\t%s\n", "Vertex", "    Dist", "    r_-", "    r_+");

  for(i=0; i<top; i++)
    {
      printf("%i\t%14.8f\t%14.8f\t%14.8f\n", i, cdist[i], 
	   knot.lenvtr[((i-1+top)%top)]/(2*tan(knot.angle[i]*.5)),
	   knot.lenvtr[i]/(2*tan(knot.angle[i]*.5)));
    }

  cout << "Smallest = " << smallestCrit << endl;
  cout << "Biggest = " << biggestCrit << endl;
}


int NewFindRate(const knots & knot, const int & i, const int & j, double & d)
{
  int flag = 0;
  int k,n;
  double s[4], m[4], a, b, vv;
  point v;
  
  v = knot.side[i] & knot.side[j];
  
  if ((vv = v * v)!= 0) {
    vv = 1/vv;
    a = (((knot[j] - knot[i]) & knot.side[j]) * v) * vv;
    b = (((knot[j] - knot[i]) & knot.side[i]) * v) * vv;
  }
  else a = b = -1;
  
  if (0 <= a && a<=1 && 0<=b && b<=1)
    {flag = 1;}
  
  else  { 
    for (n = k = 0; k < 4; k++) // find n such that m[n] is minimal
      {  m[k] = (k<2) ? min(knot, i, knot[(j+k)%knot.vnum()], s[k])
	   : min(knot, j, knot[i+k-2], s[k]);
      if (m[n] > m[k])  n = k;  
      }
    if (n<2)  a = s[n], b = n;
    else      a = n-2,  b = s[n];
  }

  d = ( (knot[i] + a*knot.side[i]) - (knot[j] + b*knot.side[j]) ).norm();

  return flag;
}

double min(const knots & knot, const int & i, const point & p, double & ratio) 
{   
  point e = p - knot[i];
  ratio = (e * knot.side[i])/(knot.side[i] * knot.side[i]);
  if (ratio < 0) ratio = 0;
  if (ratio > 1) ratio = 1;
  return (e - (ratio * knot.side[i])).norm();
}

int IsDC(const knots & knot, const int & i, const int & j, const double & a)
{
  point v=(knot[j]+a*knot.side[j]-knot[i]);

  int s;

  if(i>0)
    s=( (v * knot.side[i-1]) >= 0);
  else
    s=( (v * knot.side[knot.vnum()-1]) >= 0);
  
  int t=(v*knot.side[i] <= 0);
  if (s == t)
    return 1;
  else
    return 0;
}

void Torsion(knots & knot)
{
  double totalTorsion = 0;
  int second, third;
  cout << endl;
  int top = knot.vnum();
  vectors torsion(top);

  cout.setf(ios::fixed);
  cout.precision(6);

  totalTorsion = 0;

  for(int k=0; k<top; k++)
    {
      second = ((k+1)%top);
      third = ((k+2)%top);
      point projection(knot.side[k][0],knot.side[k][1],0);
      Normalize(knot,second,third);
      double numerator=(knot.side[k] * projection);
      double denominator = knot.side[k].norm() * projection.norm();
      double frac = numerator/denominator;

      if(frac >= 1)
	frac = 0;
      else if(frac <= -1)
	frac = M_PI;
      else
	frac = acos(frac);

      torsion[k] = frac;

      totalTorsion += torsion[k];
      cout << "torsion+[" << k << "] = " << torsion[k]  
	   << "\tt[" << k << "] = "
	   << 2*tan(torsion[k]*.5)/knot.side[k].norm() << endl;
    }
  cout << "Total torsion of knot is " << totalTorsion << endl << endl;


}

void Torsion2(knots & knot)
{
  double totalTorsion = 0;
  point cross1(0,0,0), cross2(0,0,0);
  double temp;
  cout << endl;
  int top = knot.vnum();
  vectors torsion(top);

  cout.setf(ios::fixed);
  cout.precision(6);


  int e1, e2, e3;
  totalTorsion = 0;

  for(int k=0; k<top; k++)
    {
      e1 = ((k-1+top)%top);
      e2 = k;
      e3 = ((k+1)%top);

      /* This is the torsion at edge k */
      cross1 = knot.side[e1] & knot.side[e2];
      cross2 = knot.side[e2] & knot.side[e3];

      // cout << "cross1 " << cross1 << " cross2 " << cross2 << endl;

      temp = (cross1*cross2)/(cross1.norm()*cross2.norm());

      // cout << "temp " << temp << endl;

      if(temp >= 1)
	temp = 0;
      else if(temp <= -1)
	temp = M_PI;
      else
	temp = acos(temp);

      totalTorsion += temp;

      torsion[e2] = (2*tan(temp*.5)/(knot.lenvtr[e2]));

      cout << "torsion[" << e2 << "] = "
	   << torsion[e2] << endl;
    }
  cout << "Total torsion of knot is " << totalTorsion << endl << endl;


}

void QHull(knots & knot)
{
  Resize(knot,100);
  string infile("temporaryin");
  string outfile("temporaryout");

  ofstream fout(infile.c_str());

  int top = knot.vnum();

  fout << "3\n";
  fout << top << endl;
  for(int i=0; i<top; i++){
    fout << knot[i][0] << " " << knot[i][1] << " " << knot[i][2] << endl;
  }

  fout.close();

  system("/home/rawdon/WebDownload/qhull-2003.1/src/qconvex FA < temporaryin > temporaryout");

  char line[256];
  double hsa, hv;
  string crap;

  ifstream fin(outfile.c_str());
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin.getline(line,256);
  fin >> crap >> crap >> crap >> hsa;
  fin >> crap >> crap >> hv;

  cout << "Hull Surface Area: " << hsa << endl;
  cout << "Hull Volume: " << hv << endl;
}

void ConvexHull(knots & knot)
{
  Resize(knot,100);

  int i,j,k, isbelow, isabove,l,numfaces=0;
  int itop = knot.vnum()-2, jtop = itop+1, ktop = jtop+1;
  point v, u, cross, middle(0,0,0);
  vector faces(ktop*ktop*ktop);
  double area, height,volume=0;

  for(i=0; i<itop; i++){
    for(j=i+1; j<jtop; j++){
      for(k=j+1; k<ktop; k++){
	Normalize(knot,i,j);
	point p = knot[k];
	p[0] = 0;
	double theta = acos(p[1]/ (p.norm()));
	
	if(p[2]>0)
	  theta = -theta;

	knot.rotateyz(theta,-1,ktop);

	/* v[i] is origin, v[j] is on x-axis, v[k] lies on xy plane */

	/* determine if there is a vertex lying above & one below xy plane */
	isbelow = isabove = 0;

	for(l=0; l<ktop; l++){
	  if((!isbelow) && (knot[l][2]<0)){
	    if( (l!=i) && (l!=j) && (l!=k) )
	      isbelow = 1;
	  }
	  if((!isabove) && (knot[l][2]>0)){
	    if( (l!=i) && (l!=j) && (l!=k) )
	      isabove = 1;
	  }
	}

	/* print the vertices and increment numfaces */
	if( (isbelow == 0) || (isabove == 0) )
	  {
	    faces[numfaces]=point(i,j,k);
	    numfaces ++;
	  }
      }
    }
  }

  /* determine the middle point of the knot */
  for(l=0; l<ktop; l++)
    middle = middle + knot[l];
  
  ktop = 1/ktop;
  middle[0] = middle[0]*ktop;
  middle[1] = middle[1]*ktop;
  middle[2] = middle[2]*ktop;
  
  /* sum the volumes of the pyramids determined by the face and 
     the middle point of the knot */
  for(l=0; l<numfaces; l++){
    i = int(faces[l][0]);
    j = int(faces[l][1]);
    k = int(faces[l][2]);
    
    u = (knot[j] - knot[i]);
    v = (knot[k] - knot[i]);
    
    cross = u&v;
    area = .5 * cross.norm();  /* area of base triangle */
    
    u = middle - knot[i];      /* reuse u */
    height = (fabs(u * cross) / cross.norm());
    
    /* cout << "area[" << l << "] = " << area << endl;
       cout << "height[" << l << "] = " << height << endl;
       cout << "volume[" << l << "] = " << area*height/3 << endl; */

    volume += area*height/3;
  }
  

  cout << "Total volume is " << volume << endl;
}


void Hull2(knots & knot,double & volume, double & surfacearea)
{
  Resize(knot,100);

  int i,j,k,l,isit,numfaces=0;
  int itop = knot.vnum()-2, jtop = itop+1, ktop = jtop+1;
  point v, u, cross, middle(0,0,0);
  vector faces(ktop*ktop*ktop);
  double area, height;
  volume=0;
  surfacearea=0;
  int flag = 0;

  for(i=0; i<itop; i++){
    for(j=i+1; j<jtop; j++){
      for(k=j+1; k<ktop; k++){
	cross = (knot[j]-knot[i]) & (knot[k]-knot[i]);
	l=0;
	while( (l==i) || (l==j) || (l==k) )
	  l++;
	
	isit = 1;
	double value = (knot[l]-knot[i])*cross;

	if(value<0)
	  {
	    for(l++;l<ktop; l++)
	      {
		if(((knot[l]-knot[i])*cross > 0) && (l!=i) && (l!=j) && (l!=k))
		  {
		    isit = 0;
		    break;
		  }
	      }
	  }
	else if(value>0)
	  {
	    for(l++;l<ktop; l++)
	      {
		if(((knot[l]-knot[i])*cross < 0) && (l!=i) && (l!=j) && (l!=k))
		  {
		    isit = 0;
		    break;
		  }
	      }
	  }
	else
	  flag = 1;
	
	/* print the vertices and increment numfaces */
	if(isit)
	  {
	    // cout << i << " " << j << " " << k << endl;
	    faces[numfaces]=point(i,j,k);
	    numfaces ++;
	  }
      }
    }
  }

  // cout << "middle = " << middle << endl;
  
  /* determine the middle point of the knot */
  for(l=0; l<ktop; l++){
    middle = middle + knot[l];
  }  

  double multiplier = 1/(double(ktop));

  middle[0] = middle[0]*multiplier;
  middle[1] = middle[1]*multiplier;
  middle[2] = middle[2]*multiplier;
  
  // cout << "middle = " << middle << endl;
  // cout << "numfaces = " << numfaces << endl;

  /* sum the volumes of the pyramids determined by the face and 
     the middle point of the knot */
  for(l=0; l<numfaces; l++){
    i = int(faces[l][0]);
    j = int(faces[l][1]);
    k = int(faces[l][2]);
    
    u = (knot[j] - knot[i]);
    v = (knot[k] - knot[i]);
    
    cross = u&v;
    area = .5 * cross.norm();  /* area of base triangle */
    surfacearea += area;

    u = middle - knot[i];      /* reuse u */
    height = (fabs(u * cross) / cross.norm());
    
    // cout << "middle - knot[" << i << "] = " << middle - knot[i] << endl;
    // cout << "cross = " << cross << endl;
    // cout << "height = " << height << endl;

    volume += area*height/3;

    // cout << "area = " << area << endl;
    // cout << "volume = " << volume << endl;
  }

  if(flag)
    surfacearea = volume = 0;

  cout.setf(ios::fixed);
  cout.precision(6);
  cout << "Total surface area is " << surfacearea << endl;
  cout << "Total volume is " << volume << endl;
}

void EqCriticalRun(knots & knot, const int & runs, const double & anglechange)
{
  /* this will check the two directions at a random vertex to 
     determine whether the knot is at a local min */
  
  int vnum=knot.vnum();
  int i, j, k, flag;
  knots tempknot;
  double tempRopeLength;
  int onejump=1, twojump=1;

  while( (onejump == 1) || (twojump == 1) )
    {
      onejump = twojump = 0;

      for(k=0; k<runs; k++)
	{
	  tempknot = knot;
	  tempRopeLength = knot.RopeLength();
	  
	  /* choose a random vertex */
	  i = int(drand48()*vnum);
	  j = ((i+2)%vnum);
	  flag = (i<j);  /* true as long as j doesn't wrap around */
	  
	  Normalize(knot,i,j);  /* translate and rotate */
	  
	  if(flag)
	    Crank(knot,anglechange,i,j);
	  else
	    Crank(knot,anglechange,j,i);
	  
	  knot.Compute();
	  
	  /* compareThickness(knot, tempknot) for the 1 jump at (i+1) */
	  if (knot.RopeLength() >= tempRopeLength)
	    {
	      if(flag)
		Crank(knot,-2*anglechange,i,j);
	      else
		Crank(knot,-2*anglechange,j,i);
	      
	      knot.Compute();
	      
	      if (knot.RopeLength() >= tempRopeLength)
		{
		  knot = tempknot;
		}
	      else
		{
		  onejump = 1;
		}
	    }
	  else
	    {
	      onejump = 1;
	    }
	  
	  tempknot = knot;
	  tempRopeLength = knot.RopeLength();

	  j = ((i+3)%vnum);
	  flag = (i<j);
	  
	  Normalize(knot,i,j);  /* translate and rotate */
	  
	  if(flag)
	    Crank(knot,anglechange,i,j);
	  else
	    Crank(knot,anglechange,j,i);
	  
	  knot.Compute();
	  
	  /* compareThickness(knot, tempknot) for the 1 jump at (i+1) */
	  if (knot.RopeLength() >= tempRopeLength)
	    {
	      if(flag)
	    Crank(knot,-2*anglechange,i,j);
	      else
		Crank(knot,-2*anglechange,j,i);
	      
	      knot.Compute();
	  
	      if (knot.RopeLength() >= tempRopeLength)
		{
		  knot = tempknot;
		}
	      else
		{
		  twojump = 1;
		}
	    }
	  else
	    {
	      twojump = 1;
	    }
	}
    }
}

void AutomateFindCrit(knots & knot)
{
  int runs;
  cout << "Number of iterations: ";
  cin >> runs;

  if(runs>0){
    for(double k=.1; k>=.000000000000001; k*=.1)
      EqCriticalRun(knot,runs,k);
  }
}

void NewBoxDimensions(knots & knot)
{
}
//   Resize(knot,100);

//   int i,j,k,l,m,isit,numfaces=0;
//   int itop = knot.vnum()-2, jtop = itop+1, ktop = jtop+1;
//   point v, u, cross;
//   vector faces(ktop*ktop*ktop);
//   double dtemp=0,a,bd;
//   int holder[2];
//   int segs[ktop][ktop];
//   int flag = 0;

//   for(i=0; i<ktop; i++)
//     for(j=0; j<ktop; j++)
//       segs[i][j]=0;

//   for(i=0; i<itop; i++){
//     for(j=i+1; j<jtop; j++){
//       for(k=j+1; k<ktop; k++){
// 	cross = (knot[j]-knot[i]) & (knot[k]-knot[i]);
// 	l=0;
// 	while( (l==i) || (l==j) || (l==k) )
// 	  l++;
	
// 	isit = 1;
// 	double value = (knot[l]-knot[i])*cross;

// 	if(value<0)
// 	  {
// 	    for(l++;l<ktop; l++)
// 	      {
// 		if(((knot[l]-knot[i])*cross > 0) && (l!=i) && (l!=j) && (l!=k))
// 		  {
// 		    isit = 0;
// 		    break;
// 		  }
// 	      }
// 	  }
// 	else if(value>0)
// 	  {
// 	    for(l++;l<ktop; l++)
// 	      {
// 		if(((knot[l]-knot[i])*cross < 0) && (l!=i) && (l!=j) && (l!=k))
// 		  {
// 		    isit = 0;
// 		    break;
// 		  }
// 	      }
// 	  }
// 	else
// 	  flag = 1;
	
// 	/* print the vertices and increment numfaces */
// 	if(isit)
// 	  {
// 	    faces[numfaces] = point(i,j,k);
// 	    numfaces ++;
// 	    segs[i][j] = 1;
// 	    segs[i][k] = 1;
// 	    segs[j][k] = 1;
// 	  }
//       }
//     }
//   }
  
//   flag = 0;

//   /* find maximum of minimum distances between planes and points 
//      and lines and lines */
//   for(l=0; l<numfaces; l++){
//     i = int(faces[l][0]);
//     j = int(faces[l][1]);
//     k = int(faces[l][2]);
    
//     u = (knot[j] - knot[i]);
//     v = (knot[k] - knot[i]);
    
//     cross = u&v;

//     /* planes and points */
//     for(m=0; m<ktop; m++){
//       dtemp = (fabs((knot[m]-knot[i])*cross))/cross.norm();
//       if(dtemp > d){
// 	d = dtemp;
// 	holder[0] = l;
// 	holder[1] = m;
//       }
//     }

//     /* lines and lines */
//     for(i=0; i<jtop; i++){
//       for(j=i+1; j<ktop; j++){
// 	if(segs[i][j]){
// 	  knot.findrate(i,j,a,b,dtemp);
// 	  dtemp = ((vertex[i] + a*side[i]) - (vertex[j] + b*side[j])).norm();
// 	  if(dtemp > d){
// 	    d = dtemp;
// 	    flag = 1;
// 	    holder[0] = i;
// 	    holder[1] = j;
// 	  }
// 	}
//       }
//     }

//     i = holder[0];
//     j = holder[1];
//   }

//   /* project to plane (case plane to point)*/
//   if(!flag){
    
    

//   cout.setf(ios::fixed);
//   cout.precision(6);
//   cout << "Total surface area is " << sa << endl;
//   cout << "Total volume is " << volume << endl;

//   knot.Compute();
  
//   cout << "Box dimensions: " << d1 << " X " << d2 << " X " << d3 << endl;
//   cout << "Box volume: " << d1*d2*d3 << endl;
//   cout << "Box surface area: " << 2*d1*d2 + 2*d2*d3 + 2*d1*d3 << endl;
//   cout << "sides: " << v1 << " " << v2 << endl;
//   cout << "sides: " << v3 << " " << v4 << endl;
//   cout << "sides: " << v5 << " " << v6 << endl;
// }


void IsItCrit(knots & knot)
{  
/* this will check the two directions at all of the vertices to 
     determine whether the knot is at a local min */
  
  int vnum=knot.vnum();
  int i, j, flag;
  knots tempknot;
  int yesorno = 1;
  double tempRopeLength;
  double anglechange = 0.000000000000001;

  cout.setf(ios::fixed);
  cout.precision(16);

  tempknot = knot;
  tempRopeLength = knot.RopeLength();
	  
  for(i=0; i<vnum; i++){
    j = ((i+2)%vnum);
    flag = (i<j);  /* true as long as j doesn't wrap around */
	  
    Normalize(knot,i,j);  /* translate and rotate */
    
    if(flag)
      Crank(knot,anglechange,i,j);
    else
      Crank(knot,anglechange,j,i);
    
    knot.Compute();
    
    /* compareThickness(knot, tempknot) for the 1 jump at (i+1) */
    if (knot.RopeLength() < tempRopeLength){
      cout << "original = " << tempRopeLength << endl;
      cout << "new = " << knot.RopeLength() << endl;
      cout << "difference = " << tempRopeLength - knot.RopeLength() << endl;
      yesorno = 0;
      cout << "i = " << i << "  j = " << j << endl;
    }

    knot = tempknot;
    
    j = ((i+3)%vnum);
    flag = (i<j);
	  
    Normalize(knot,i,j);  /* translate and rotate */
    
    if(flag)
      Crank(knot,anglechange,i,j);
    else
      Crank(knot,anglechange,j,i);
	  
    knot.Compute();
	  
    /* compareThickness(knot, tempknot) for the 2 jump at (i+1) */
    if (knot.RopeLength() < tempRopeLength){
      cout << "original = " << tempRopeLength << endl;
      cout << "new = " << knot.RopeLength() << endl;
      cout << "difference = " << tempRopeLength - knot.RopeLength() << endl;
      yesorno = 0;
      cout << "i = " << i << "  j = " << j << endl;
    }
  }

  knot = tempknot;

if (yesorno == 1)
  cout << "Sure as shootin', it is critical.\n";
else
  cout << "Sorry, this is not a critical point.  Try again.\n";
}

void CritFind(knots & knot)
{  
/* this will check the two directions at all of the vertices to 
     determine whether the knot is at a local min */
  
  int vnum=knot.vnum();
  int i, j, flag;
  int counter;
  knots tempknot;
  int yesorno = 0;

  double tempRopeLength;
  double anglechange = 0.1;

  cout.setf(ios::fixed);
  cout.precision(16);

  tempknot = knot;
  tempRopeLength = knot.RopeLength();
	  
  for(;anglechange >= 0.000000000000001; anglechange*=0.1){
    counter = 0;
    while(yesorno==0){
      counter++;
      yesorno = 1;
      for(i=0; i<vnum; i++){
	j = ((i+2)%vnum);
	flag = (i<j);  /* true as long as j doesn't wrap around */
	
	Normalize(knot,i,j);  /* translate and rotate */
	
	if(flag)
	  Crank(knot,anglechange,i,j);
	else
	  Crank(knot,anglechange,j,i);
	
	knot.Compute();
	
	/* compareThickness(knot, tempknot) for the 1 jump at (i+1) */
	if (knot.RopeLength() < tempRopeLength){
	  yesorno = 0;
	  tempRopeLength = knot.RopeLength();
	  tempknot = knot;
	}
	else{
	  if(flag)
	    Crank(knot,-2*anglechange,i,j);
	  else
	    Crank(knot,-2*anglechange,j,i);
	  
	  knot.Compute();
	  
	  /* compareThickness(knot, tempknot) for the 1 jump at (i+1) */
	  if (knot.RopeLength() < tempRopeLength){
	    yesorno = 0;
	    tempRopeLength = knot.RopeLength();
	    tempknot = knot;
	  }
	  else
	    knot = tempknot;
	}

	j = ((i+3)%vnum);
	flag = (i<j);
	
	Normalize(knot,i,j);  /* translate and rotate */
	
	if(flag)
	  Crank(knot,anglechange,i,j);
	else
	  Crank(knot,anglechange,j,i);
	
	knot.Compute();
	
	/* compareThickness(knot, tempknot) for the 2 jump at (i+1) */
	if (knot.RopeLength() < tempRopeLength){
	  yesorno = 0;
	  tempRopeLength = knot.RopeLength();
	  tempknot = knot;
	}
	else{
	  if(flag)
	    Crank(knot,-2*anglechange,i,j);
	  else
	    Crank(knot,-2*anglechange,j,i);
	  
	  knot.Compute();
	  
	  /* compareThickness(knot, tempknot) for the 1 jump at (i+1) */
	  if (knot.RopeLength() < tempRopeLength){
	    yesorno = 0;
	    tempRopeLength = knot.RopeLength();
	    tempknot = knot;
	  }
	  else
	    knot = tempknot;
	}
      }
    }   
    Save(knot,"mostrecent");
    yesorno = 0;
  } 
}

void NonEqCritFind(knots & knot)
{  
/* this will check the three directions at all of the vertices to 
     determine whether the knot is at a local min. 
   This is for non-equilateral knots. */
  
  int vnum=knot.vnum();
  int i, j, flag;
  int counter;
  knots tempknot;
  int yesorno = 0;

  double tempRopeLength;
  double change = 0.1;

  cout.setf(ios::fixed);
  cout.precision(16);

  tempknot = knot;
  tempRopeLength = knot.RopeLength();
	  
  cout << "Here we go\n";

  for(;change >= 0.000000000000001; change*=0.1){
    counter = 0;
    while(yesorno==0){
      counter++;
      yesorno = 1;
      for(i=0; i<vnum; i++){
	for(j=0; j<3; j++){
	  
	  knot[i][j] = knot[i][j] + change;
	  knot.Compute();
	
	  /* compareThickness(knot, tempknot) for the + jump */
	  if (knot.RopeLength() < tempRopeLength){
	    yesorno = 0;
	    tempRopeLength = knot.RopeLength();
	    tempknot = knot;
	  }
	  else{
	    knot[i][j] = knot[i][j] - 2*change;
	    knot.Compute();
	    
	    /* compareThickness(knot, tempknot) for the - jump */
	    if (knot.RopeLength() < tempRopeLength){
	      yesorno = 0;
	      tempRopeLength = knot.RopeLength();
	      tempknot = knot;
	    }
	    else
	      knot = tempknot;
	  }
	}
      }
    }
    Save(knot,"mostrecent");
    yesorno = 0;
  } 
}

void Distance(const knots & knot)
{
  cout << "Load what other knot: ";

  knots knot2;
  string filename;

  cin >> filename;
  
  ifstream fin(filename.c_str());

  fin >> knot2;

  fin.close();

  if(knot.vnum() != knot2.vnum()){
    cout << "Can't compare knots with different numbers of edge.\n";
    return;
  }

  double sum=0;
  point p;

  for(int i=0; i<knot.vnum(); i++){
    p = knot[i] - knot2[i];
    cout << p << endl;
    sum += p.norm();
  }

  cout << "Distance = " << sum << endl;
}

void IsItEq(knots & knot)
{
  Resize(knot, knot.vnum());

  double mu = knot.injrad();
  
  double shortestside = knot.Length();
  double longestside = 0;

  for(int i=0; i<knot.vnum(); i++){
    if(knot.lenvtr[i] < shortestside)
      shortestside = knot.lenvtr[i];
    if(knot.lenvtr[i] > longestside)
      longestside = knot.lenvtr[i];
  }

  double spread = longestside - 1;
  if(1 - shortestside < spread)
    spread = 1 - shortestside;

  if(shortestside < mu)
    mu = shortestside;

  double q1 = mu/knot.vnum();
  double q2 = (mu*mu)/4;

  cout << "mu/n = " << q1 << endl;
  cout << "mu^2/4 = " << q2 << endl;
  cout << "| L - 1 | = " << spread << endl;

  if((spread<q1) && (spread<q2))
    cout << "These computations show that an equilateral exists!\n";
  else
    cout << "Sorry, we cannot conclude anything.\n";
}
  
/* translates knot so that vertex_0 is at the origin and vertex_1 lies
   on the positive x-axis and v_3 lies in first quadrant*/
void DefaultPosition(knots & knot)
{
  int k, top = knot.vnum();

  Normalize(knot,0,1);

  cout << knot << endl;

  double angle = knot[2][1] / ( sqrt( knot[2][1]*knot[2][1] +
                                           knot[2][2]*knot[2][2]));

  if(angle >= 1)
    angle = 0;
  else if(angle <= -1)
    angle = M_PI;
  else
    angle = acos(angle);

  cout << knot[2] << endl;
  cout << angle << endl;

  if(knot[2][2] >= 0)
    angle = angle;
  else
    angle = -angle;

  knot.rotateyz(angle);

  cout << knot << endl;
}

void Mobius(knots & knot)
{
  cout << "The transformation is of the form (ax+b)/(cx+d).  "
       << "Enter a, b, c, and d: ";
  double a,b,c,d;
  cin >> a >> b >> c >> d;

  cout << "Tranformation is (" << a << "x + " << b << ") / ("
       << c << "x + " << d << ")\n";

  int top=knot.vnum(), i, j;
  double coord;

  if((a*d-b*c) != 0){
    for(i=0; i<top; i++){
      for(j=0; j<3; j++){
	coord = knot[i][j];
	knot[i][j] = (a*coord+b)/(c*coord+d);
      }
    }
  }
  else
    cout << "Bad determinant!\n";
}
    
void IsDCincharge(const knots & knot)
{
  double minrad, ii, vi, vv, dcsd;
  minrad = knot.MinRad();
  ii = knot.InteriorInterior();
  vi = knot.VertexInterior();
  vv = knot.VertexVertex();

  dcsd = ii;
  if(vi < dcsd)
    dcsd = vi;
  if(vv < dcsd)
    dcsd = vv;
  
  if(minrad < dcsd)
    cout << "MinRad is too small.\n";
  else
    cout << "dcsd is in charge, hurrah!\n";
}


/* this is a lie right now */
/* I am confident this works (except for planar knots) on May 15, 2003 */
/* Except I haven't implemented looking for line line pairs for which
   we could get the height.  So I am assuming in general that plane point
   pairs bound the minimum height. */
void SkinnyMaxBoxDimensions(knots & knot, double & length, double & width, double & height)
{
  /* find the exterior vertices of the convex hull */
  Resize(knot,100);

  int i,j,k,l,isit;
  int itop = knot.vnum()-2, jtop = itop+1, ktop = jtop+1;
  const int vnumsquared = ktop*ktop;
  point v, u, cross, middle(0,0,0);
  // vectors isitonoutside(ktop,0);
  int flag;

  vector outsideplanes(itop*itop*itop);
  int planecounter = 0;

  int lengthvertex;
  length = 0.0;
  int lengthi, lengthj, lengthk;


  for(i=0; i<itop; i++){
    for(j=i+1; j<jtop; j++){
      for(k=j+1; k<ktop; k++){
	flag = 0;
	
	cross = (knot[j]-knot[i]) & (knot[k]-knot[i]);
	l=0;
	while( (l==i) || (l==j) || (l==k) )
	  l++;
	
	isit = 1;
	double value = (knot[l]-knot[i])*cross;

	if(value<0)
	  {
	    for(l++;l<ktop; l++)
	      {
		if(((knot[l]-knot[i])*cross > 0) && (l!=i) && (l!=j) && (l!=k))
		  {
		    isit = 0;
		    break;
		  }
	      }
	  }
	else if(value>0)
	  {
	    for(l++;l<ktop; l++)
	      {
		if(((knot[l]-knot[i])*cross < 0) && (l!=i) && (l!=j) && (l!=k))
		  {
		    isit = 0;
		    break;
		  }
	      }
	  }
	else
	  flag = 1;
	
	/* print the vertices */
	if(isit)
	  {
	    outsideplanes[planecounter] = point(i,j,k);
	    planecounter++;
	    /* cout << i << " " << j << " " << k << endl; */
	  }
      }
    }
  }

  /* for(k=0; k<ktop; k++){
     cout << k << " " << isitonoutside[k] << endl;
     } */

  /* search for the slimmest planes enclosing the knot */

  /* first look for a three point plane versus one point situation */
  double projdist;
  double smallest = 650000;
  

  /* if the three bound a plane on the exterior
     then search for point that maxes distance to plane */
  double tempmaxdist = 0.0;
  int tempmaxdistvertex;
  int m;

  for(m=0; m<planecounter; m++){
 
    tempmaxdist = 0.0;

    i = int(outsideplanes[m][0]);
    j = int(outsideplanes[m][1]);
    k = int(outsideplanes[m][2]);

 
    /* cout << "i " << i << " j " << j << " k " << k << endl; */
    
    for(int n=0; n<ktop; n++){
      if( (n!= i) && (n!= j) && (n!=k) ){
	cross = ((knot[j] - knot[i]) & (knot[k] - knot[i])); 
	projdist = (fabs((cross  * (knot[n] - knot[i]))))/cross.norm();
	
	if(projdist > tempmaxdist){
	  /* cout << "i " << i << " j " << j << " k " << k << " n " << n 
	     << " okay " << endl; 
	     cout << "tempmaxdist " << tempmaxdist << " projdist " << projdist
	     << endl;
	  */
	  
	  tempmaxdist = projdist;
	  tempmaxdistvertex = n;
	}
	else if(projdist == tempmaxdist){
	  cerr << "problems in computing length!  Two pairs with same distance\n";
	  cerr << "i " << i << " j " << j << " k " << k << " n " << n << endl;
	  exit(0);
	}
      }
    }
    
    /* compare to old data to find largest distance between
       bounding 3 point planes versus single points */
    
    if(tempmaxdist > length){
      length = tempmaxdist;
      lengthvertex = tempmaxdistvertex;
      lengthi = i;
      lengthj = j;
      lengthk = k;
      
      /* cout << "length " << length << " vertex " << lengthvertex <<
	 " i " << lengthi << " j " << lengthj << " k " << lengthk << endl; */
    }
  }

  cout << "Final length " << length << " vertex " << lengthvertex
       << " i " << lengthi << " j " << lengthj << " k " << lengthk << endl;


  /* Note Dec 16, 2003 */
  /* I think the line line search is going to be fruitless */
  /* It seems non-generic and time consuming, so I say screw it */
  /* Now look for line/line pairs maximizing distance */
  /* int islineonoutside[ktop][ktop];

  for(int i=0; i<itop; i++)
  islineonoutside = 0;

  for(int i=0; i<planecounter; i++){
    islineonoutside[outsideplanes[i][0]][outsideplanes[i][1]] = 1;
    islineonoutside[outsideplanes[i][0]][outsideplanes[i][2]] = 1;
    islineonoutside[outsideplanes[i][1]][outsideplanes[i][2]] = 1;
  }

  double maxlines = 0.0;
  double linedist;

  for(int j=0; j<jtop; j++){
    for(int k=j+1; k<ktop; k++){
      if(islineonoutside[j][k] == 1){
	for(int l=j; l<jtop; l++){
	  for(int m=0; m<ktop; m++){
	    if(islineonoutside[l][m] == 1){
	      minifindrate(i,j,l,m,linedist);
	      if(linedist > maxlines)
		maxlines = linedist;
	    }
	  }
	}
      }
    }
  }
	    
  */

  /* project onto the plane normal to the triple point plane by moving the
     lengthvertex to origin making the plane parallel to the yz-plane */
  cross = ((knot[lengthj] - knot[lengthi]) & (knot[lengthk] - knot[lengthi])); 
  PointsNormalize(knot,knot[lengthvertex],knot[lengthvertex]-cross);
  
  knots projection(knot);

  for(i=0; i<ktop; i++)
    projection[i][0] = 0;

  /* check to see that this is working
  cross = ((knot[lengthj] - knot[lengthi]) & (knot[lengthk] - knot[lengthi])); 
  cout << "cross is " << cross << endl;
  cout << ((projection[lengthj] - projection[lengthi]) & 
  (projection[lengthk] - projection[lengthi])) << endl;
  it appears to be okay */

  /* Save(knot,"knot");
     Save(projection,"projknot"); */

  /* find fattest line/point pair, that is, the set of lines bounding the
     projection that are closest together */ 

  /* if all of the cross products of (vk-vi) with (vj-vi) are in the same
     direction, then the line between vi and vj is on the outside of
     the convex hull of the planar figure */

  point tempcross;
  width=0.0;
  int widthi, widthj, widthvertex;
  tempmaxdist = 0.0;
  // int tempmaxdistvertex, m;

  for(i=0; i<jtop; i++){ /* start of line */
    for(j=i+1; j<ktop; j++){ /* end of line */
      k=0;

      while( (k==i) || (k==j) )
	k++;
      cross = (projection[k]-projection[i]) & (projection[j] - projection[i]);

      for(k++; k<ktop; k++){ /* other point */
	isit = 1;
	if( (k!=i) && (k!=j) ){
	  tempcross = (projection[k]-projection[i]) & (projection[j] - projection[i]);
	  if((tempcross * cross) < 0){
	    isit = 0;
	    /* cout << "i " << i << " j " << j << " k " << k << " bad\n"; */
	    break;
	  }
	}
      }

      point linedirection;

      if(isit == 1){
	tempmaxdist = 0.0;
	/* cout << "i " << i << " j " << j << " good\n"; */

	for(m=0; m<ktop; m++){
	  if( (m!=i) && (m!=j) && (m!=k) ){
	    linedirection = (projection[i] - projection[j]);
	    projdist = ((projection[m] - projection[i])&linedirection).norm()/linedirection.norm();
	    if(tempmaxdist < projdist){
	      tempmaxdist = projdist;
	      tempmaxdistvertex = m;
	    }
	    else if(projdist == tempmaxdist){
	      cerr << "problems in computing width!\n";
	      exit(0);
	    }
	  }
	}

	/* the width is the smallest of these */

	if(tempmaxdist > width){
	  width = tempmaxdist;
	  widthvertex = tempmaxdistvertex;
	  widthi = i;
	  widthj = j;
	  
	  /* cout << "width " << width << " vertex " << widthvertex
	     << " widthi " << widthi << " widthj " << widthj << endl; */
	}
      }
    }
  }

  cout << "Final width " << width << " vertex " << widthvertex
       << " widthi " << widthi << " widthj " << widthj << endl;

  /* project again, this time on the perpendicular to chord vivj */
  
  cross = projection[widthj] - projection[widthi];
  /* cout << "vector is " << cross << endl; */
  cross[0] = 0;

  double theta = acos(cross[1]/(cross.norm()));
  if(cross[2]<0){
    theta = -theta;
  }

  projection.rotateyz(theta);
  /* Save(projection,"zprojknot"); */

  /* project onto y-axis */

  int heighti, heightj;
  double miny=650000.0, maxy=-650000.0;

  for(i=0; i<ktop; i++){
    projection[i][2] = 0;
    if(projection[i][1] < miny){
      miny = projection[i][1];
      heighti = i;
    }
    if(projection[i][1] > maxy){
      maxy = projection[i][1];
      heightj = i;
    }
  }
  
  height = maxy - miny;

  cout << "Final height " << height << " i " << heighti 
       << " j " << heightj << endl;
  
  /* Save(projection,"zzprojknot"); */

  cout << "Max Box Dimensions: length " << length << " width " << width 
       << " height " << height << endl;
  cout << "Max Box Surface Area: " 
       << 2*length*width+2*length*height+2*width*height << endl;
  cout << "Max Box Volume: " << length*width*height << endl;

  if(length < width)
    cout << "Oops: in max box, length < width\n";
  if(width < height)
    cout << "Oops: in max box, width < height\n";

}


/* I am confident this works (except for planar knots) on May 15, 2003 */
/* Except I haven't implemented looking for line line pairs for which
   we could get the height.  So I am assuming in general that plane point
   pairs bound the minimum height. */
void SkinnyMinBoxDimensions(knots & knot, double & length, double & width, double & height)
{
  /* find the exterior vertices of the convex hull */
  Resize(knot,100);

  int i,j,k,l,isit;
  int itop = knot.vnum()-2, jtop = itop+1, ktop = jtop+1;
  point v, u, cross, middle(0,0,0);
  // vectors isitonoutside(ktop,0);
  int flag;

  vector outsideplanes(itop*itop*itop);
  int planecounter = 0;

  int heightvertex;
  height = 65000.0;
  int heighti, heightj, heightk;


  for(i=0; i<itop; i++){
    for(j=i+1; j<jtop; j++){
      for(k=j+1; k<ktop; k++){
	flag = 0;
	
	cross = (knot[j]-knot[i]) & (knot[k]-knot[i]);
	l=0;
	while( (l==i) || (l==j) || (l==k) )
	  l++;
	
	isit = 1;
	double value = (knot[l]-knot[i])*cross;

	if(value<0)
	  {
	    for(l++;l<ktop; l++)
	      {
		if(((knot[l]-knot[i])*cross > 0) && (l!=i) && (l!=j) && (l!=k))
		  {
		    isit = 0;
		    break;
		  }
	      }
	  }
	else if(value>0)
	  {
	    for(l++;l<ktop; l++)
	      {
		if(((knot[l]-knot[i])*cross < 0) && (l!=i) && (l!=j) && (l!=k))
		  {
		    isit = 0;
		    break;
		  }
	      }
	  }
	else
	  flag = 1;
	
	/* print the vertices */
	if(isit)
	  {
	    outsideplanes[planecounter] = point(i,j,k);
	    planecounter++;
	    /* cout << i << " " << j << " " << k << endl; */
	  }
      }
    }
  }

  /* for(k=0; k<ktop; k++){
     cout << k << " " << isitonoutside[k] << endl;
     } */

  /* search for the slimmest planes enclosing the knot */

  /* first look for a three point plane versus one point situation */
  double projdist;
  double smallest = 650000;
  

  /* if the three bound a plane on the exterior
     then search for point that maxes distance to plane */
  double tempmaxdist = 0.0;
  int tempmaxdistvertex;
  int m;

  for(m=0; m<planecounter; m++){
 
    tempmaxdist = 0.0;

    i = int(outsideplanes[m][0]);
    j = int(outsideplanes[m][1]);
    k = int(outsideplanes[m][2]);

 
    /* cout << "i " << i << " j " << j << " k " << k << endl; */
    
    for(int n=0; n<ktop; n++){
      if( (n!= i) && (n!= j) && (n!=k) ){
	cross = ((knot[j] - knot[i]) & (knot[k] - knot[i])); 
	projdist = (fabs((cross  * (knot[n] - knot[i]))))/cross.norm();
	
	if(projdist > tempmaxdist){
	  /* cout << "i " << i << " j " << j << " k " << k << " n " << n 
	       << " okay " << endl; 
	  cout << "tempmaxdist " << tempmaxdist << " projdist " << projdist
	  << endl; */
	  
	  tempmaxdist = projdist;
	  tempmaxdistvertex = n;
	}
	else if(projdist == tempmaxdist){
	  cerr << "problems in computing height!\n";
	  cerr << "i " << i << " j " << j << " k " << k << " n " << n << endl;
	  exit(0);
	}
      }
    }
    
    /* compare to old data to find smallest distance between
       bounding 3 point planes versus single points */
    
    if(tempmaxdist < height){
      height = tempmaxdist;
      heightvertex = tempmaxdistvertex;
      heighti = i;
      heightj = j;
      heightk = k;
      
      /* cout << "height " << height << " vertex " << heightvertex <<
	 " i " << heighti << " j " << heightj << " k " << heightk << endl; */
    }
  }

  cout << "Final height " << height << " vertex " << heightvertex
       << " i " << heighti << " j " << heightj << " k " << heightk << endl;


  /* Now looking for line-line pairs that might bound the darn thing */









  /* project onto the plane normal to the triple point plane by moving the
     heightvertex to origin making the plane parallel to the yz-plane */
  cross = ((knot[heightj] - knot[heighti]) & (knot[heightk] - knot[heighti])); 
  PointsNormalize(knot,knot[heightvertex],knot[heightvertex]-cross);
  
  knots projection(knot);

  for(i=0; i<ktop; i++)
    projection[i][0] = 0;

  /* check to see that this is working
  cross = ((knot[heightj] - knot[heighti]) & (knot[heightk] - knot[heighti])); 
  cout << "cross is " << cross << endl;
  cout << ((projection[heightj] - projection[heighti]) & 
  (projection[heightk] - projection[heighti])) << endl;
  it appears to be okay */

  /* Save(knot,"knot");
     Save(projection,"projknot"); */

  /* find slimmest line/point pair, that is, the set of lines bounding the
     projection that are closest together */ 

  /* if all of the cross products of (vk-vi) with (vj-vi) are in the same
     direction, then the line between vi and vj is on the outside of
     the convex hull of the planar figure */

  point tempcross;
  width=650000.0;
  int widthi, widthj, widthvertex;
  tempmaxdist = 0.0;
  // int tempmaxdistvertex, m;

  for(i=0; i<jtop; i++){ /* start of line */
    for(j=i+1; j<ktop; j++){ /* end of line */
      k=0;

      while( (k==i) || (k==j) )
	k++;
      cross = (projection[k]-projection[i]) & (projection[j] - projection[i]);

      for(k++; k<ktop; k++){ /* other point */
	isit = 1;
	if( (k!=i) && (k!=j) ){
	  tempcross = (projection[k]-projection[i]) & (projection[j] - projection[i]);
	  if((tempcross * cross) < 0){
	    isit = 0;
	    /* cout << "i " << i << " j " << j << " k " << k << " bad\n"; */
	    break;
	  }
	}
      }

      point linedirection;

      if(isit == 1){
	tempmaxdist = 0.0;
	/* cout << "i " << i << " j " << j << " good\n"; */

	for(m=0; m<ktop; m++){
	  if( (m!=i) && (m!=j) && (m!=k) ){
	    linedirection = (projection[i] - projection[j]);
	    projdist = ((projection[m] - projection[i])&linedirection).norm()/linedirection.norm();
	    if(tempmaxdist < projdist){
	      tempmaxdist = projdist;
	      tempmaxdistvertex = m;
	    }
	    else if(projdist == tempmaxdist){
	      cerr << "problems in computing width!\n";
	      exit(0);
	    }
	  }
	}

	/* the width is the smallest of these */

	if(tempmaxdist < width){
	  width = tempmaxdist;
	  widthvertex = tempmaxdistvertex;
	  widthi = i;
	  widthj = j;
	  
	  /* cout << "width " << width << " vertex " << widthvertex
	     << " widthi " << widthi << " widthj " << widthj << endl; */
	}
      }
    }
  }

  cout << "Final width " << width << " vertex " << widthvertex
       << " widthi " << widthi << " widthj " << widthj << endl;

  /* project again, this time on the perpendicular to chord vivj */
  
  cross = projection[widthj] - projection[widthi];
  /* cout << "vector is " << cross << endl; */
  cross[0] = 0;

  double theta = acos(cross[1]/(cross.norm()));
  if(cross[2]<0){
    theta = -theta;
  }

  projection.rotateyz(theta);
  /* Save(projection,"zprojknot"); */

  /* project onto y-axis */

  int lengthi, lengthj;
  double miny=650000.0, maxy=-650000.0;

  for(i=0; i<ktop; i++){
    projection[i][2] = 0;
    if(projection[i][1] < miny){
      miny = projection[i][1];
      lengthi = i;
    }
    if(projection[i][1] > maxy){
      maxy = projection[i][1];
      lengthj = i;
    }
  }
  
  length = maxy - miny;

  cout << "Final length " << length << " i " << lengthi 
       << " j " << lengthj << endl;
  
  /* Save(projection,"zzprojknot"); */

  cout << "Min Box Dimensions: length " << length << " width " << width 
       << " height " << height << endl;
  cout << "Min Box Surface Area: " 
       << 2*length*width+2*length*height+2*width*height << endl;
  cout << "Min Box Volume: " << length*width*height << endl;

  if(length < width)
    cout << "Oops: in min box, length < width\n";
  if(width < height)
    cout << "Oops: in min box, width < height\n";
}

 
void Report(knots & knot, const string & file)
{
  double a, b, c;

  /* Convex Hull Data */
  Hull2(knot,a,b);
  cout << "Convex Hull\n";
  cout << "Hull Surface Area = " << b << endl;
  cout << "Hull Volume = " << a << endl;

  /* Standard Box Data */
  Standardize(knot,a,b,c);
  cout << "Standard Box\n";
  cout << "Box Length = " << a << endl;
  cout << "Box Width = " << b << endl;
  cout << "Box Height = " << c << endl;
  cout << "Box Surface Area = " << 2*(a*b+a*c+b*c) << endl;
  cout << "Box Volume = " << a*b*c << endl;
  
  /* MinBox Data */
  SkinnyMinBoxDimensions(knot,a,b,c);
  cout << "Min Box\n";
  cout << "Box Length = " << a << endl;
  cout << "Box Width = " << b << endl;
  cout << "Box Height = " << c << endl;
  cout << "Box Surface Area = " << 2*(a*b+a*c+b*c) << endl;
  cout << "Box Volume = " << a*b*c << endl;

  /* MinBox Data */
  SkinnyMaxBoxDimensions(knot,a,b,c);
  cout << "Max Box\n";
  cout << "Box Length = " << a << endl;
  cout << "Box Width = " << b << endl;
  cout << "Box Height = " << c << endl;
  cout << "Box Surface Area = " << 2*(a*b+a*c+b*c) << endl;
  cout << "Box Volume = " << a*b*c << endl;

  /* Ropelength Data */
  knot.Display();

 
}

void RadiusOfGyration(const knots & knot)
{
  /* middle is the center of mass */
  point middle(0.0,0.0,0.0);
  int top = knot.vnum();

  for(int i=0; i<top; i++){
    middle = middle + knot[i];
  }

  /* We're going to need this a couple of times, so instead of doing
     the type casting over and over, I'm just going to hit it once */
  double invtop = 1.0/double(top);

  /* okay, now have the center of mass */
  middle = invtop*middle;

  /* radgyr is the Mean Square Radius of Gyration, this means take the 
     squared distance to the center of mass and average it */
  double radgyr = 0;

  /* this is here for historic reasons, I love history */
  double tmp;

  for(int i=0; i<top; i++){
    /* distance squared is the dot product with itself */
    tmp = (knot[i]-middle)*(knot[i]-middle);
    radgyr += tmp;
  }

  /* Take the average of the squared distances */
  radgyr = radgyr*invtop;

  /* this is all one word so that other scripts I have written still work */
  cout << "MeanSquareRadius of gyration = " << radgyr << endl;
}


int MiniBall(const knots & knot)
{
  const int numpoints = knot.vnum();
  Miniball<3> mb;

  Miniball<3>::Point p;
  for(int i=0; i<numpoints; i++){
    p = knot[i];
    mb.check_in(p);
  }

  // construct ball, using the pivoting method
  // -----------------------------------------
  cout << "Constructing miniball..."; cout.flush();
  mb.build();
  cout << "done." << endl << endl;
  
  // output center and squared radius
  // --------------------------------
  cout << "Center:         " << mb.center() << endl;
  cout << "Squared radius: " << mb.squared_radius() << endl << endl;
  
  // output number of support points
  // -------------------------------
  cout << mb.nr_support_points() << " support points: " << endl << endl;
   
  // output support points
  // ---------------------
  Miniball<3>::Cit it;
  for (it=mb.support_points_begin(); it!=mb.support_points_end(); ++it)
    cout << *it << endl;
  cout << endl;
  
  // output accuracy
  // ---------------
  double slack;
  cout << "Relative accuracy: " << mb.accuracy (slack) << endl;
  cout << "Optimality slack:  " << slack << endl;
  
  return 0;
}


void CurvaturePlot(knots & knot){
  string file;
  cout << "File for curvature output: ";
  cin >> file;

  ofstream fout(file.c_str());

  int top = knot.vnum();

  for(int k=0; k<top; k++){
    fout << k << "\t" << knot.lenvtr[((k-1+top)%top)]/(2*tan(knot.angle[k]*.5))/knot.injrad()
	 << "\t" << knot.lenvtr[k]/(2*tan(knot.angle[k]*.5))/knot.injrad() << endl;
  }
  
  fout.close();
}

void MakeThicknessOne(knots & knot){
  double ratio = 1.0/knot.injrad();

  int top = knot.vnum();

  for(int i=0; i<top; i++)
    knot[i] = ratio*knot[i];
  
  knot.Compute();
}

double Asphericity(const knots & knot){

  double matrixentries[9];
  double multiplisor;
  double asphericity = 0.0;
  double sum, tmp;
  double evalm[3];
  int edges = knot.vnum();

  // http://www.gnu.org/software/gsl/manual/gsl-ref_14.html
  gsl_matrix_view mat = gsl_matrix_view_array(matrixentries,3,3);
  gsl_vector *eval = gsl_vector_alloc(3);
  
  // Can comment out if don't care about axes
  gsl_matrix *evec = gsl_matrix_alloc(3,3);

  // gsl_eigen_symm_workspace * ws = gsl_eigen_symm_alloc(3);
  gsl_eigen_symmv_workspace * ws = gsl_eigen_symmv_alloc(3);
  

  multiplisor = 1.0/(2.0*edges*edges);
  asphericity = 0.0;
  
  for(int row=0; row<3; row++){
    for(int col=row; col<3; col++){
      
      sum = 0.0;
      
      for(int i=0; i<edges; i++){
	for(int j=0; j<edges; j++){
	  sum += (knot[i][row]-knot[j][row])*(knot[i][col]-knot[j][col]);
	}
      }
      
      tmp = multiplisor * sum;
      matrixentries[3*col+row] = tmp;
      matrixentries[3*row+col] = tmp;
    }
  }
  
  /* use symmv instead of symm to be able to get eigenvectors */
  // gsl_eigen_symm(&mat.matrix, eval, ws);
  gsl_eigen_symmv(&mat.matrix, eval, evec, ws);

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  for(int m=0; m<3; m++){
    evalm[m] = gsl_vector_get(eval,m);
    gsl_vector_view evec_i = gsl_matrix_column(evec,m);

    printf("eigenvalue = %f\n",evalm[m]);
    printf("eigenvector = \n");
    gsl_vector_fprintf(stdout,&evec_i.vector,"%g");
  }
  
  asphericity += 0.5*(pow(evalm[0]-evalm[1],2) + pow(evalm[1]-evalm[2],2) + pow(evalm[2]-evalm[0],2))/(pow(evalm[0]+evalm[1]+evalm[2],2)); 
  
  gsl_eigen_symmv_free(ws);

  return asphericity;

}

double MDEnergy(const knots & knot){


}


double SymmetryScore(const knots & knot){

  double score, group_avg, dthisvert, vstep, score_contrib;
  int nverts = knot.vnum();
  int nvertsover2 = nverts/2;
  double dnverts = double(nverts);
  int thisvert, v, w;

  double top = 20.0;
  
  if(dnverts < top)
    top = dnverts/2;;

  for(double rottype=2; rottype<top; rottype++){ // type of symmetry, test up to 20
    double members[int(rottype)];
    score = 0;

    for(int s=1; s<=nvertsover2; s++) { // Skip is set to s.
      
      vstep = dnverts/rottype;
      
      for(int i=0; i<nverts; i++) {
	
	group_avg = 0;
	
	thisvert = i;
	dthisvert = double(i);
	for(int j=0; j<rottype; j++) {
	  
	  v = thisvert % nverts;
	  w = (thisvert + s) % nverts;
	  
	  members[j] = (knot[v] - knot[w]).norm();
	  
	  group_avg += members[j];
	  dthisvert += vstep;
	  thisvert = int(dthisvert);
	}

	group_avg /= rottype;  // to get the average
	
	score_contrib=0;
	for(int j=0; j<rottype; j++) {
	  score_contrib += fabs(members[j] - group_avg)/(double)(rottype);
	}
	
	score += score_contrib;
      }
    }

    cout << int(rottype) << "-fold symmetry score: " << score << endl;
  }

  return 1.0;
}


