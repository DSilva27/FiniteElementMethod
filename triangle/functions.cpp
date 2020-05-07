#include <iostream>
#include <vector>
#include "triangle.h"

using namespace std;


Triangle::Triangle(vector <vector <double>> vert,
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

