#include "LInt.h"

/*
 Here we apply the fact that the variables of a line integral can be parametrized in terms of a parameter t.

 For a straight line joining the points (a,b) and (a',b'). The parametrization is as simple as:

 x = (a' - a)*t + a; y = (b' - b)*t + b; for 0 <= t <= 1.

 The integrand has the general expression:

 f(x(t), y(t)) * [ (a' - a)^2 + (b' - b)^2 ]^0.5, where the second terms is the Jacobian of the transformation.

 */
LineIntegrator::LineIntegrator(void){}

//Calculates the Jacobian for two given points
double LineIntegrator::dJ( double p1[2], double p2[2]){

  return sqrt( pow(( p2[0] - p1[0] ),2) + pow(( p2[1] - p2[1] ),2) );
}

//Calculates x or y in terms of t
double LineIntegrator::Transf( double t, double c1, double c2 ) {return t*( c2 - c1 ) + c1;}

//Transforms the integrand to the parametrized form and calculates it
double LineIntegrator::TransfIntegrand( double t, double (*integrand)(double, double), double p1[2], double p2[2] ){

  double x = Transf( t, p1[0], p2[0]);
  double y = Transf( t, p1[1], p2[1]);

  double I = (*integrand)( x, y );
  I *= dJ( p1, p2 );

  return I;
}

//n must be a positive, even integer
double LineIntegrator::LineIntegral( double (*integrand)(double,double),
                                     double p1[2], double p2[2], double n ){

  //Here we integrate using the Composite Simpson's rule
  double h = 1/n;

  double XI0 = TransfIntegrand(0, (*integrand), p1, p2) + TransfIntegrand(1, (*integrand), p1, p2);
  double XI1 = 0;
  double XI2 = 0;

  for (int i = 1; i < n; i++){

    double X = i*h;

    if (i%2 == 0){

      XI2 += TransfIntegrand(X, (*integrand), p1, p2);
    }

    else {
      XI1 += TransfIntegrand(X, (*integrand), p1, p2);
    }
  }

  double XI = h*( XI0 + 2*XI2 + 4*XI1)/3;

  return XI;
}
