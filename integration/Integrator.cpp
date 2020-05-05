#include "Integrator.h"

TriangleIntegrator::TriangleIntegrator(void){}

double TriangleIntegrator::dJ( double u, double v, double p1[2], double p2[2], double p3[2]){

  double dxdu = ( (1 - v) * p2[0] + v * p2[0] - p1[0] );
  double dxdv = ( u*p3[0] - u*p2[0] );

  double dydu = ( (1 - v) * p2[1] + v * p2[1] - p1[1] );
  double dydv = ( u*p3[1] - u*p2[1] );

  return std::abs( dxdu*dydv - dxdv*dydu );
}

double TriangleIntegrator::Transf( double u, double v, double c1, double c2, double c3 ) {return (1 - u)*c1 + u*( (1 - v)*c2 + v*c3 );}

double TriangleIntegrator::TransfIntegrand( double (*integrand)(double, double), double u, double v, double p1[2], double p2[2], double p3[2] ){

  double x = Transf( u, v, p1[0], p2[0], p3[0]);
  double y = Transf( u, v, p1[1], p2[1], p3[1]);

  double I = (*integrand)( x, y );
  I *= dJ( u, v, p1, p2, p3 );

  return I;
}

double TriangleIntegrator::DoubleIntegral(double (*integrand)(double,double),
                                          double p1[2], double p2[2], double p3[2],
                                          double h, double k
                                          ){

  double lx = 0;
  double ux = 1;
  double ly = 0;
  double uy = 1;

  int nx, ny;
  
  // z stores the table
  // ax[] stores the integral wrt y
  // for all x points considered
  double z[500][500];
  double ax[500];

  double answer;

  // Calculating the numner of points
  // in x and y integral
  nx = (ux - lx) / h + 1;
  ny = (uy - ly) / k + 1;

  // Calculating the values of the table
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      z[i][j] = TransfIntegrand( (*integrand), lx + i * h, ly + j * k, p1, p2, p3);
    }
  }

  // Calculating the integral value
  // wrt y at each point for x
  for (int i = 0; i < nx; ++i) {
    ax[i] = 0;
    for (int j = 0; j < ny; ++j) {

      if (j == 0 || j == ny - 1)
        ax[i] += z[i][j];
      else if (j % 2 == 0)
        ax[i] += 2 * z[i][j];
      else
        ax[i] += 4 * z[i][j];
    }
    ax[i] *= (k / 3);
  }

  answer = 0;

  // Calculating the final integral value
  // using the integral obtained in the above step
  for (int i = 0; i < nx; ++i) {
    if (i == 0 || i == nx - 1)
      answer += ax[i];
    else if (i % 2 == 0)
      answer += 2 * ax[i];
    else
      answer += 4 * ax[i];
  }
  answer *= (h / 3);

  return answer;
}
