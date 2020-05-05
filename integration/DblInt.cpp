//Here we use the fact that an integral over any triangle can be transformed easily into an integral over a square (0,1)x(0,1)
//When we are integrating over the triangle the variables are x,y
//When we are integrating over the square the variables are u, v
#include "DblInt.h"

TriangleIntegrator::TriangleIntegrator(void){}

//Calculates the Jacobian of the transformation
double TriangleIntegrator::dJ( double u, double v, double p1[2], double p2[2], double p3[2]){

  double dxdu = ( (1 - v) * p2[0] + v * p2[0] - p1[0] );
  double dxdv = ( u*p3[0] - u*p2[0] );

  double dydu = ( (1 - v) * p2[1] + v * p2[1] - p1[1] );
  double dydv = ( u*p3[1] - u*p2[1] );

  return std::abs( dxdu*dydv - dxdv*dydu );
}

//Calculates x or y, given u and v
double TriangleIntegrator::Transf( double u, double v, double c1, double c2, double c3 ) {return (1 - u)*c1 + u*( (1 - v)*c2 + v*c3 );}

//Calculates the integrand for a given u and v
double TriangleIntegrator::TransfIntegrand( double (*integrand)(double, double), double u, double v, double p1[2], double p2[2], double p3[2] ){

  double x = Transf( u, v, p1[0], p2[0], p3[0]);
  double y = Transf( u, v, p1[1], p2[1], p3[1]);

  double I = (*integrand)( x, y );
  I *= dJ( u, v, p1, p2, p3 );

  return I;
}

//Calculates the double integral using Simpson's Rule
double TriangleIntegrator::DoubleIntegral(double (*integrand)(double,double), //Integrand
                                          double p1[2], double p2[2], double p3[2], //Vertices
                                          double n, double m //Parameters for integration
                                          ){

  //Here we apply Simpson's Double Integral method. This is not a general method. It is simplified for the current problem
  //in order to save memory and optimize the algorithm

  double h = 1./n;
  double J1 = 0.;
  double J2 = 0.;
  double J3 = 0.;

  for( int i = 0; i <= n; i++ ){

    double x = i*h;
    double HX = 1./m;

    double K1 = TransfIntegrand( (*integrand), x, 0., p1, p2, p3) + TransfIntegrand( (*integrand), x, 1., p1, p2, p3);
    double K2 = 0.;
    double K3 = 0.;


    for (int j = 1; j < m; j++){

      double y = j*HX;
      double Q = TransfIntegrand( (*integrand), x, y, p1, p2, p3);

      if (j%2 == 0) K2 += Q;

      else K3 += Q;
    }

    double L = (K1 + 2.*K2 + 4.*K3)*HX/3.;

    if (i == 0 || i == n) J1 += L;

    else if( i%2 == 0) J2 += L;

    else J3 += L;
  }

  double J = h*( J1 + 2.*J2 + 4.*J3 ) / 3.;

  return J;
}
