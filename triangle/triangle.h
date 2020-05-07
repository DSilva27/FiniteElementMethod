#ifndef TRIANGLE_ELEMENT_H
#define TRIANGLE_ELEMENT_H

#include <iostream>
#include <vector>

using namespace std;

class Triangle{
    friend class FiniteElement;
    
    public:
    Triangle(vector <vector <double>>, vector <int>, vector <int>);
    ~Triangle();
    
    //private:
    vector <vector <double>> vertices;
    vector <int> nodes;
    vector <int> boundary;
};
