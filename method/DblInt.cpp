//Here we use the fact that an integral over any triangle can be transformed easily into an integral over a square (0,1)x(0,1)
//When we are integrating over the triangle the variables are x,y
//When we are integrating over the square the variables are u, v
#include "Math/DblInt.h"

TriangleIntegrator::TriangleIntegrator(void){}

//Calculates the Jacobian of the transformation
double TriangleIntegrator::dJ( double u, double v, mat vert ){

  double dxdu = ( (1 - v) * vert[1][0] + v * vert[1][0] - vert[0][0] );
  double dxdv = ( u*vert[2][0] - u*vert[1][0] );

  double dydu = ( (1 - v) * vert[1][1] + v * vert[1][1] - vert[0][1] );
  double dydv = ( u*vert[2][1] - u*vert[1][1] );

  return std::abs( dxdu*dydv - dxdv*dydu );
}

//Calculates x or y, given u and v
double TriangleIntegrator::Transf( double u, double v, double c1, double c2, double c3 ) {return (1 - u)*c1 + u*( (1 - v)*c2 + v*c3 );}

//Calculates the integrand for a given u and v
double TriangleIntegrator::TransfIntegrand( func f, mat vert, double u, double v ){

  double x = Transf( u, v, vert[0][0], vert[1][0], vert[2][0]);
  double y = Transf( u, v, vert[0][1], vert[1][1], vert[2][1]);

  double I = f(x,y);

  I *= dJ( u, v, vert );

  return I;
}

double TriangleIntegrator::TransfIntegrand( func f, vec N1, mat vert, double u, double v ){

  double x = Transf( u, v, vert[0][0], vert[1][0], vert[2][0]);
  double y = Transf( u, v, vert[0][1], vert[1][1], vert[2][1]);

  double I;

  I = f(x, y) * ( N1[0] + N1[1]*x + N1[2]*y );

  I *= dJ( u, v, vert );

  return I;
}

double TriangleIntegrator::TransfIntegrand( func f, vec N1, vec N2, mat vert, double u, double v ){

  double x = Transf( u, v, vert[0][0], vert[1][0], vert[2][0]);
  double y = Transf( u, v, vert[0][1], vert[1][1], vert[2][1]);

  double I;

  I = f(x, y) * ( N1[0] + N1[1]*x + N1[2]*y ) * ( N2[0] + N2[1]*x + N2[2]*y );

  I *= dJ( u, v, vert );

  return I;
}

//Calculates the double integral using Simpson's Rule
double TriangleIntegrator::DoubleIntegral(func f, //Integrand
                                          mat vert, //Vertices
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

    double K1 = TransfIntegrand( f, vert, x, 0. ) + TransfIntegrand( f, vert, x, 1. );
    double K2 = 0.;
    double K3 = 0.;


    for (int j = 1; j < m; j++){

      double y = j*HX;
      double Q = TransfIntegrand( f, vert, x, y );

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

double TriangleIntegrator::DoubleIntegral(func f, //Integrand
                                          vec N1,
                                          mat vert, //Vertices
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

    double K1 = TransfIntegrand( f, N1, vert, x, 0. ) + TransfIntegrand( f, N1, vert, x, 1. );
    double K2 = 0.;
    double K3 = 0.;


    for (int j = 1; j < m; j++){

      double y = j*HX;
      double Q = TransfIntegrand( f, N1, vert, x, y );

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

double TriangleIntegrator::DoubleIntegral(func f, //Integrand
                                          vec N1,
                                          vec N2,
                                          mat vert, //Vertices
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

    double K1 = TransfIntegrand( f, N1, N2, vert, x, 0. ) + TransfIntegrand( f, N1, N2, vert, x, 1. );
    double K2 = 0.;
    double K3 = 0.;


    for (int j = 1; j < m; j++){

      double y = j*HX;
      double Q = TransfIntegrand( f, N1, N2, vert, x, y );

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
