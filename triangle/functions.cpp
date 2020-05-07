#include <iostream>
#include <vector>
#include "triangle.h"

using namespace std;


Triangle::Triangle(vector <vector <double>> vert,
                   vector <int> nod){
    cout << "Bienvenido" << endl;
    
    vertices = vert;
    nodes = nod;
}


Triangle::~Triangle(){
    cout << "Adios" << endl;
}

