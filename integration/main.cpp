#include "Integrator.h"
#include <iostream>

double Integrand( double, double );

int main(){

  double P1[2] = {3.,1.};
  double P2[2] = {2.,2.};
  double P3[2] = {4.,2.};

  TriangleIntegrator TI;;

  cout << TI.DoubleIntegral(Integrand, P1, P2, P3, 0.005, 0.005) << endl;

  return 0;
}

double Integrand(double x, double y){
  return x*y;
}

