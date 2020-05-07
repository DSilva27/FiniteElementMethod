#include <iostream>
#include <vector>
#include "finite_element.h"

using namespace std;

int main(){
  
  FiniteElement ex;
  
  ex.load_data();
  
  for (int i=0; i<ex.M; i++){
    
    cout << ex.elements[i].boundary[2] << endl;
    
  }
  
  
  return 0;
}

