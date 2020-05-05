#include <iostream>
#include <vector>

using namespace std;

class FiniteElement{
    
    public:
    FiniteElement();
    void solve();
    ~FiniteElement();
    
    private:
    int K;
    int N;
    int M;
    int n;
    int m;
    vector <vector <double>> vertex;
};
