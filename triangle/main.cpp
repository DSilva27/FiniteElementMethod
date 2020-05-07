#include <iostream>
#include <vector>
#include "triangle.h"

using namespace std;

int main(){
    vector <vector <double>> vec1(3, vector <double> (2));
    vector <vector <double>> vec2(3, vector <double> (2));
    vector <vector <double>> vec3(3, vector <double> (2));
    vector <int> nod1(3, 1), nod2(3, 6), nod3(3, 355);
    
    for (int i=0; i<3; i++){
        vec1[i][0] = 0.1;
        vec1[i][1] = 0.4567;
        
        vec2[i][0] = 654.1;
        vec2[i][1] = 277.4567;
        
        vec3[i][0] = 2.76;
        vec3[i][1] = 54.7;
    }
    
    /*
    Triangle ex1(vec1, nod1), ex2(vec2, nod2), ex3(vec3, nod3);
    
    for (int i=0; i<3; i++){
        cout << ex2.vertices[i][0] << endl;
        cout << ex2.nodes[i] << endl;
    }
    */
    vector < Triangle > pr;
    vector <vector <double>> vec;
    vector <int> nod, bound;
    vector <int> nodp(3, 10);
    
    
    for (int i=0; i<3; i++){
        
        //vec.push_back({1, 2});
        vec.push_back({1, 2});
        nod.push_back(4);
        bound.push_back(1);
        
    }
    
    for (int i=0; i<3; i++){
        if (i==1){ 
            Triangle ex(vec, nod, bound);
            pr.push_back(ex);
        }
        
        else{
            Triangle ex(vec, nodp, bound);
            pr.push_back(ex);
        }
    }
    
    
    cout << pr[0].nodes[1] << endl;
    cout << pr[1].nodes[1] << endl;
    
    
    
    
    
    return 0;
}

