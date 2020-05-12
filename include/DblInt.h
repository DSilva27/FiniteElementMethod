#ifndef DBLINT_H
#define DBLINT_H

#include <cmath>
#include <vector>

using namespace std;

typedef double ( *func )( double, double );
typedef vector< double > vec;
typedef vector< vec > mat;

class TriangleIntegrator{

 public:
     
  TriangleIntegrator( void );

  //Calculates the Jacobian for the transformation from a triangle to a square
  double dJ( double, double, mat );

  //Transforms from the square coordinates to the original cartesian coordinates
  double Transf( double, double, double, double, double );

  //Evaluates the transformed integrand

  double TransfIntegrand( func, mat, double, double );

  double TransfIntegrand( func, vec, mat, 
                          double, double );

  double TransfIntegrand( func, vec, vec, mat,
                          double, double );


  //Calculates the double integral from certain integrand in a triangle defined by three vertices

  double DoubleIntegral( func, mat, int, int );

  double DoubleIntegral( func, vec, mat, 
                         int, int );

  double DoubleIntegral( func, vec, vec, mat, int, int );
};

#endif
