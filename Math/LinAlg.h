#ifndef LINALG_H
#define LINALG_H

#include <iostream>
//#include <armadillo>
#include <vector>
#include <cmath>

using namespace std;
//using namespace arma;

typedef vector < double > vec;
typedef vector < vec > mat;

class Linalg{
 public:
  Linalg(void);

  //  void Solve( vector< double > &, vector< vector< double > >, vector< double > );

  void SOR(const mat, const vec, vec &, const double, const double, const int);
  double DistBetVectors( const vec, const vec);

  //double Determinant( vector< vector< double > > );
  double Det33( mat );
};


#endif
