#include <iostream>
#include <vector>
#include "finite_element.h"

using namespace std;

//Functions appearing at the left side
double P(double, double);
double Q(double, double);
double R(double, double);

double F(double, double);

double G(double, double);
double G1(double, double);
double G2(double, double);

/*
  uxx + uyy = 0
  g(x,y) = 0
 */

int main(){

  vfunc VFUNC{ P, Q, R, F, G, G1, G2 };
  FiniteElement ex;
  vector< double > gamma;
  vector< vector< vector < double > > > N_coef;
  double res;
  
  ex.load_data();
  
  ex.solve( VFUNC );
  
  ex.results_to_variable( gamma, N_coef );
  
  ex.generate_data();
  
  return 0;
}


double P(double x, double y){ return 1; }

double Q(double x, double y){ return 1; }

double R(double x, double y){ return 0; }

double F(double x, double y){ return 0; }

double G(double x, double y){ return 4; }

double G1(double x, double y){ return 0; }

double G2(double x, double y){

  if ( ((x >= 0 && x <= 0.2) && (y >= 0.2 && y <= 0.4)) or ((x > 0.4 && x <= 0.5) && (y >= 0.1 && y <= 0.2)) ){
    
      return (x + y)/2;
  }
  
  else if ( ((x > 0.2 && x <= 0.4) && (y == 0.2)) or ((x > 0.5 && x <= 0.6) && (y == 0.1 )) ){
  
    return x;
  }
  
  else if( x == 0.6 && (y >= 0 && y < 0.1)){
  
    return y;
  }
  
  return 0;
}
