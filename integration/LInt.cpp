#include "LInt.h"

LineIntegrator::LineIntegrator(void){}

double LineIntegrator::dJ( double p1[2], double p2[2]){

  return sqrt( pow(( p2[0] - p1[0] ),2) + pow(( p2[1] - p2[1] ),2) );
}

double LineIntegrator::Transf( double t, double c1, double c2 ) {return t*( c2 - c1 ) + c1;}

double LineIntegrator::TransfIntegrand( double t, double (*integrand)(double, double), double p1[2], double p2[2] ){

  double x = Transf( t, p1[0], p2[0]);
  double y = Transf( t, p1[1], p2[1]);

  double I = (*integrand)( x, y );
  I *= dJ( p1, p2 );

  return I;
}

double LineIntegrator::LineIntegral( double (*integrand)(double,double),
                                          double p1[2], double p2[2],
                                          double n ){

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
