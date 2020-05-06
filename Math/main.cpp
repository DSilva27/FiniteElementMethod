#include "DblInt.h"
#include "LInt.h"
#include <iostream>

double Integrand( double, double );
double Integrand2( double, double );
double Integrand3( double, double );
int main(){

  double P1[2] = {3.,1.};
  double P2[2] = {2.,2.};
  double P3[2] = {4.,2.};

  vector <double (*)(double, double)> v { Integrand, Integrand2 };
  TriangleIntegrator TI;
  LineIntegrator LT;

  cout << TI.DoubleIntegral(v, P1, P2, P3, 1000, 1000) << endl;
  cout << TI.DoubleIntegral(Integrand3, P1, P2, P3, 1000, 1000) << endl;
  cout << LT.LineIntegral(v, P2, P3, 1000) << endl;;
  return 0;
}

double Integrand(double x, double y){ return y; }

double Integrand2(double x, double y){ return x; }

double Integrand3(double x, double y){ return x*y; }
