#ifndef LINALG_H
#define LINALG_H

#include <iostream>
//#include <armadillo>
#include <vector>
#include <cmath>

using namespace std;
//using namespace arma;


class Linalg{
 public:
  Linalg(void);

  //  void Solve( vector< double > &, vector< vector< double > >, vector< double > );

  void SOR(const vector < vector< double > >, const vector < double >, vector < double > &, const double, const double, const int);
  double DistBetVectors( vector < double >, vector < double >);

  //double Determinant( vector< vector< double > > );
  double Det33(vector <vector <double> >);
};


#endif
