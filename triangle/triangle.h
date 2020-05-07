#include <iostream>
#include <vector>

using namespace std;

class Triangle{
    friend class FiniteElement;
    
    public:
    Triangle(vector <vector <double>>, vector <int>);
    ~Triangle();
    
    private:
    vector <vector <double>> vertices;
    vector <int> nodes;
};
