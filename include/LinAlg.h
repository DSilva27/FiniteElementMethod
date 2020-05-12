#ifndef LINALG_H
#define LINALG_H

#include <iostream>

#include <vector>
#include <cmath>

using namespace std;

typedef vector< double > vec;
typedef vector< vec > mat;

class Linalg{
    
 public:
     
  Linalg( void );

  void SOR( mat, vec, vec &,  double,
            double,  int );

  double DistBetVectors(  vec,  vec );

  double Det33( mat );
};


#endif
