#include <iostream>
#include <vector>
#include "triangle.h"

using namespace std;


Triangle::Triangle(mat vert,
                   vector <int> nod){
  cout << "Bienvenido" << endl;
  
  vertices = vert;
  nodes = nod;
}


Triangle::~Triangle(){
  cout << "Adios" << endl;
}

