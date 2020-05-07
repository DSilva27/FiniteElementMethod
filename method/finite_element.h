#include <iostream>
#include <vector>

#include "triangle.h"

using namespace std;

typedef vector< double > vec;
typedef vector< vec > mat;
typedef vector< mat > cube;

class FiniteElement{
  
  public:
  FiniteElement();
  void load_data();
  void solve();
  ~FiniteElement();
  
  //private:
  int K;
  int N;
  int M;
  int n;
  int m;
  vector <class Triangle> elements; // test this
};
