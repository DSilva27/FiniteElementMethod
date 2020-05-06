#ifndef DBLINT_H
#define DBLLINT_H

#include <cmath>
#include <vector>

using namespace std;

class TriangleIntegrator{

 public:
  TriangleIntegrator();

  //Calculates the Jacobian for the transformation from a triangle to a square
  double dJ( double, double, double *, double *, double *);

  //Transforms from the square coordinates to the original cartesian coordinates
  double Transf( double, double, double, double, double );

  //Evaluates the transformed integrand
  double TransfIntegrand(vector <double (*)(double, double)>, double, double, double *, double *, double * );

  //Calculates the double integral from certain integrand in a triangle defined by three vertices
  double DoubleIntegral( vector <double (*)(double, double)>, double *, double *, double *, double, double );
};

#endif
