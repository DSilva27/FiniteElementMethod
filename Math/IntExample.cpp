#include "DblInt.h"
#include "LInt.h"
#include <iostream>

double Integrand( double, double );
double Integrand2( double, double );
double Integrand3( double, double );


int main(){

  mat Vert{ {3.,1.}, {2.,2.}, {4.,2.}};
  vec N1{ 3., 3., 3.};

  TriangleIntegrator TI;
  LineIntegrator LT;

  cout << TI.DoubleIntegral(Integrand3, Vert, 1000, 1000) << endl;
  // cout << TI.DoubleIntegral(Integrand3, P1, P2, P3, 1000, 1000) << endl;
  cout << LT.LineIntegral(Integrand3, N1, Vert[0], Vert[1], 1000) << endl;
  return 0;
}

double Integrand(double x, double y){ return y; }

double Integrand2(double x, double y){ return x; }

double Integrand3(double x, double y){ return x*y; }
