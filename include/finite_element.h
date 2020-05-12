#include <iostream>
#include <vector>

#include "triangle.h"
#include "DblInt.h"
#include "LInt.h"
#include "LinAlg.h"

using namespace std;

typedef double (*func)( double, double );
typedef vector< func > vfunc;

typedef vector< double > vec;
typedef vector< vec > mat;
typedef vector< mat > cube;

/*
 Each element of vfunc will be:

 0: p(x,y)
 1: q(x,y)
 2: r(x,y)
 3: f(x,y)
 4: g(x,y)
 5: g1(x,y)
 6: g2(x,y)
 */

class FiniteElement{
  
  public:
  FiniteElement();
  void load_data();
  void solve( vfunc );
  void gamma_to_txt();
  void N_coef_to_txt();
  void results_to_txt();
  void gamma_to_variable( vec& );
  void N_coef_to_variable( cube& );
  void results_to_variable( vec&, cube& );
  ~FiniteElement();
  
  private:
  int K;
  int N;
  int M;
  int n;
  int m;
  int p;
  vector < class Triangle > elements;
  mat nodes;
  cube N_coef;
  vec gamma;

  //Math tools  
  TriangleIntegrator DInt;
  LineIntegrator LInt;
  Linalg LinAlg;
};
