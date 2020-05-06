#ifndef LINALG_H
#define LINALG_H

#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;


class Linalg{
 public:
  Linalg(void);

  void Solve( vector< double > &, vector< vector< double > >, vector< double > );

  double Determinant( vector< vector< double > > );
};


#endif
