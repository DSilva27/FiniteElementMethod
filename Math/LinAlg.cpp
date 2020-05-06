#include "LinAlg.h"

Linalg::Linalg(void){}

void Linalg::Solve( vector < double > &x, vector< vector< double > > m, vector< double > b ){

  mat A(m.size(),m.size());

  vec B(b);

  for (int i = 0; i < m.size(); i++){
    for (int j = 0; j< m.size(); j++){
      A.at(i,j) = m[i][j];
    }
  }

  vec X = solve(A,B);

  for (int k = 0; k < m.size(); k++){
    x[k] = X.at(k);
  }
}

double Linalg::Determinant(vector< vector< double > > m){

  mat A(m.size(),m.size());

  for (int i = 0; i < m.size(); i++){
    for (int j = 0; j< m.size(); j++){
      A.at(i,j) = m[i][j];
    }
  }

  return det(A);
}
