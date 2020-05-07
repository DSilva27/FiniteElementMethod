#ifndef LINT_H
#define LINT_H

#include <cmath>
#include <vector>

using namespace std;

typedef double (*func)(double, double);
typedef vector< func > vfunc;

class LineIntegrator{

 public:
  LineIntegrator();

  //Calculates the Jacobian for the transformation from a variables x,y to parameter t
  double dJ( double *, double * );

  //Calculates x(t) and y(t). The parameters are: t, c1, c2. Where c is either x or y
  double Transf( double, double, double );

  //Evaluates the integrand given t and the endpoints
  double TransfIntegrand(double, vfunc, double *, double *);
  //Calculates the line integral using the Composite Simpson's Rule
  double LineIntegral( vfunc, double *, double *, double );

};

#endif
