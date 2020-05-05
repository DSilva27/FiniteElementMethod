#include "INPUT/data_reader.h"

void ExtractInfo(int *KPTR,
		 int *NPTR,
		 int *MPTR,
		 int *nPTR,
		 int *mPTR,
		 vector<vector <double>> &vec){

  std::ifstream File;
  File.open("INPUT/triangles.txt");

  //For reading the info in File
  int node;
  double x;
  double y;

  //Declaring INPUT values for the current problem
  
  *KPTR = 2;  //Two interior triangles
  *NPTR = 6;  //Four triangles with at least an edge in S2
  *MPTR = 10; //Four leftover triangles
  *nPTR = 5;  //Five nodes in D U S2
  *mPTR = 11; //Six nodes on S1

  //Initializes the vector that will be filled
  for (int i = 0; i < *mPTR; i++){

    vector<double> subVector(2);
    vec.push_back(subVector);
  }
  
  while (File >> node >> x >> y){

    //Fills the vertex matrix
    vec[node][0] = x;
    vec[node][1] = y;
  }
}
