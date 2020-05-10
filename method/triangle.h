#ifndef TRIANGLE_ELEMENT_H
#define TRIANGLE_ELEMENT_H

#include <iostream>
#include <vector>

using namespace std;

class Triangle;
typedef vector< vector< double >> mat;

class Triangle{
  friend class FiniteElement;
  
  public:
  Triangle(mat, vector <int>);
  ~Triangle();
  
  private:
  mat vertices;
  vector <int> nodes;
};

#endif
