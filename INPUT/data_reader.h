#ifndef TRIANGLE_READER_H
#define TRIANGLE_READER_H

//#define FILE="triangles.txt"

#include <fstream>
#include <string>
#include <vector>

using namespace std;

typedef vector< double > vec;
typedef vector< vec > mat;
//This function extracts the information from the triangles file
//It assigns the corresponding x, y values for certain node
void ExtractInfo(int*, //K
        int*, //N
        int*, //M
        int*, //n
        int*, //m
        mat & //vertex
        );
  
#endif
