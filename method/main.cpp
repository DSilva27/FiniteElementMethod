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
  
  ex.load_data();
  
  /*
  for (int i=0; i<ex.M; i++){
    
    cout << ex.nodes[i][2] << endl;
    //cout << ex.elements[i].boundary[2] << endl;
    
  }*/
 
  return 0;
}


double P(double x, double y){ return 1; }

double Q(double x, double y){ return 1; }

double R(double x, double y){ return 0; }

double F(double x, double y){ return 0; }

double G(double x, double y){ return 4; }

double G1(double x, double y){ return 4; }
double G2(double x, double y){ return 4; }
