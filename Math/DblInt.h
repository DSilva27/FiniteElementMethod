#ifndef DBLINT_H
#define DBLLINT_H

#include <cmath>
#include <vector>

using namespace std;

typedef double (*func)(double, double);
typedef vector< func > vfunc;

class TriangleIntegrator{

 public:
  TriangleIntegrator();

  //Calculates the Jacobian for the transformation from a triangle to a square
  double dJ( double, double, double *, double *, double *);

  //Transforms from the square coordinates to the original cartesian coordinates
  double Transf( double, double, double, double, double );

  //Evaluates the transformed integrand
  double TransfIntegrand(vfunc, double, double, double *, double *, double * );

  //If there is just one function
  double TransfIntegrand(func, double, double, double *, double *, double * );

  //Calculates the double integral from certain integrand in a triangle defined by three vertices
  double DoubleIntegral( vfunc, double *, double *, double *, double, double );

  double DoubleIntegral( func, double *, double *, double *, double, double );
};

#endif
