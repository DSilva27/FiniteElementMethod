#include "DblInt.h"
#include "LInt.h"
#include <iostream>

double Integrand( double, double );

int main(){

  double P1[2] = {3.,1.};
  double P2[2] = {2.,2.};
  double P3[2] = {4.,2.};

  TriangleIntegrator TI;
  LineIntegrator LT;

  cout << TI.DoubleIntegral(Integrand, P1, P2, P3, 1000, 1000) << endl;
  cout << LT.LineIntegral(Integrand, P2, P3, 1000) << endl;;
  return 0;
}

double Integrand(double x, double y){
  return x*y;
}

