#include "LinAlg.h"

int main(){

  vector <vector <double>> M {{1,2},{3,4}};

  vector <double> B {5,6};

  vector <double> x(B.size(), 0);

  Linalg Solver;

  Solver.Solve(x, M, B);

  for(int i = 0; i < B.size(); i++){
    cout << x[i] << endl;
  }

  cout << Solver.Determinant(M) << endl;
  return 0;
}
