#ifndef DBLINT_H
#define DBLINT_H

#include <cmath>
#include <vector>
#include <gs1/gs1_integration.h>
#include <stdio.h>
using namespace std;

class LineIntegrator{

 public:
  LineIntegrator();

  //Calculates the Jacobian for the transformation from a triangle to a square
  double dJ( double *, double * );

  //Transforms from the square coordinates to the original cartesian coordinates
  double Transf( double, double, double );

  //Evaluates the transformed integrand
  double TransfIntegrand(double, double (*)(double, double), double *, double *);

  //Calculates the double integral from certain integrand in a triangle defined by three vertices
  double LineIntegral( double (*)(double, double), double *, double *, double );
};

#endif
