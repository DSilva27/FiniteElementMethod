#include "LinAlg.h"

Linalg::Linalg(void){}

// void Linalg::Solve( vector < double > &x, vector< vector< double > > m, vector< double > b ){

//   mat A(m.size(),m.size());

//   vec B(b);

//   for (int i = 0; i < m.size(); i++){
//     for (int j = 0; j< m.size(); j++){
//       A.at(i,j) = m[i][j];
//     }
//   }

//   vec X = solve(A,B);

//   for (int k = 0; k < m.size(); k++){
//     x[k] = X.at(k);
//   }
// }

void Linalg::SOR( const mat M, const vec b, vec &XO,
         const double w, const double TOL, const int N ){

  int k = 1;
  int n = M.size();
  vector< double > x(n, 0.0);

  while (k <= N){

    for (int i = 1; i <= n; i++){

      double sum1 = 0;

      for (int j = 1; j <= i - 1; j++){
        sum1 += M[i-1][j-1]*x[j-1];
      }

      double sum2 = 0;

      for (int j = i+1; j <= n; j++){

        sum2 += M[i-1][j-1]*XO[j-1];
      }

      x[i-1] = (1-w)*XO[i-1] + (1/M[i-1][i-1]) * ( w * ( -sum1 - sum2 + b[i-1] ) );
    }

    if (DistBetVectors( x, XO ) < TOL){

      XO = x;

      return; //To exit the function
    }

    k++;


    for (int i = 1; i <= n; i++){

      XO[i-1] = x[i-1];
    }

  }

  XO = x;

  cout << "Max number of iterations exceeded" << endl;
  return;
}

double Linalg::DistBetVectors( vec v1, vec v2){

  double dist;
  for(int i = 0; i < v1.size(); i++){

    dist += pow(v1[i] - v2[i], 2);
  }

  //  cout << dist << endl;

  return sqrt(dist);
}


// double Linalg::Determinant(vector< vector< double > > m){

//   mat A(m.size(),m.size());

//   for (int i = 0; i < m.size(); i++){
//     for (int j = 0; j< m.size(); j++){
//       A.at(i,j) = m[i][j];
//     }
//   }

//   return det(A);
// }

double Linalg::Det33(mat M){

  double det;

  det = ( M[1][1]*M[2][2] - M[2][1]*M[1][2] ) - M[0][1]*( M[2][2] - M[1][2] ) + M[0][2]*( M[2][1] - M[1][1] );

  return det;
}
