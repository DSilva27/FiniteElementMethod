#ifndef LINT_H
#define LINT_H

#include <cmath>
#include <vector>

using namespace std;

typedef double (*func)(double, double);
typedef vector< double > vec;

class LineIntegrator{

 public:
  LineIntegrator();

  //Calculates the Jacobian for the transformation from a variables x,y to parameter t
  double dJ( vec, vec );

  //Calculates x(t) and y(t). The parameters are: t, c1, c2. Where c is either x or y
  double Transf( double, double, double );

  //Evaluates the integrand given t and the endpoints
  double TransfIntegrand(func, vec, vec, vec, double);

  double TransfIntegrand(func, vec, vec, vec, vec, double);

  //Calculates the line integral using the Composite Simpson's Rule
  double LineIntegral( func, vec, vec, vec, int);

  double LineIntegral( func, vec, vec, vec, vec, int);

};

#endif
