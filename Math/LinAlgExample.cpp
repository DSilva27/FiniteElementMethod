#include "LinAlg.h"

int main(){

  vector <vector <double>> M {{1, -1, 1},{1, 6, 2}, {1, 3, 7}};

  vector <double> B {1, 0, 4};

  vector <double> x(B.size(), 0);

  Linalg Solver;

  Solver.SOR(M, B, x, 1.25, 0.003, 20);

  for(int i = 0; i < B.size(); i++){
    cout << x[i] << endl;
  }

  cout << Solver.Det33(M) << endl;
  return 0;
}
