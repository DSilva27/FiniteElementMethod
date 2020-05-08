#include <iostream>
#include <vector>
#include "triangle.h"

using namespace std;


Triangle::Triangle(mat vert,
                   vector <int> nod,
                   vector <int> bound){
  cout << "Bienvenido" << endl;
  
  vertices = vert;
  nodes = nod;
  boundary = bound;
}


Triangle::~Triangle(){
  cout << "Adios" << endl;
}

