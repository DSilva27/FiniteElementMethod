#include <iostream>
#include <vector>
#include <fstream>
#include "triangle.h" //??
#include "finite_element.h"

using namespace std;

FiniteElement::FiniteElement(){
  cout << "Bienvenido" << endl;
}


FiniteElement::~FiniteElement(){
  cout << "Adios" << endl;
}


void FiniteElement::load_data(){
  
  vector <vector <double>> vertices_vec;
  vector <int> nodes_vec;
  vector <int> boundary_vec;
  
  int count = 1; // To know in which triangle we are at
  int number, node, bound;
  double x, y;
  
  K = 2;
  N = 6;
  M = 10;
  n = 5;
  m = 11;
  
  file.open("data/data_triangles.txt");
  
  while (file >> number >> node >> x >> y << bound){
    if ( number == count ){
      
      // Fill vectors with data
      vertices_vec.push_back({x, y});
      nodes_vec.push_back(node);
      boundary_vec.push_back(boundary);
    }
    
    else{
      
      // When no more data about that triangle, push back element and increase count
      Triangle triangle(vertices_vec, nodes_vec, boundary_vec);
      elements.push_back(triangle);
      
      vertices_vec.clear();
      nodes_vec.clear();
      boundary_vec.clear();
      
      count += 1;
    }
  }
  
  file.close();
}


void FiniteElement::solve(){
  
  vector <double> gamma(m);
  vector <double> beta(n);
  vector <vector <double>> alpha(n, vector <double> (n));
  vector <vector <vector <double>>> N_coef(M, vector < vector <double>> (3, vector <double> (3)));
  vector <vector <vector <double>>> z(M, vector < vector <double>> (3, vector <double> (3)));
  vector <vector <vector <double>>> J((N-K-1), vector < vector <double>> (3, vector <double> (3)));
  vector <vector <double>> I((N-K-1), vector <double> (3));
  vector <vector <double>> H(M, vector <double> (3));
  
  double det;
  double integral_p, integral_q, integral_r;
  
  double x1, x2, x3;
  double y1, y2, y3;
  
  int l, t;
  
  
  // Step 1
  for (int l=n+1; l<m; l++){
    gamma[l] = g(vertex[l][0], vertex[l][0]); // g def is missing
  }
  
  // Step 2 is not necessary because vectors are already initialized to 0
  
  // Step 3
  for (int i=0; i<M; i++){
    det = 1; // det function is missing
    
    // Coefficients of the function N(x, y)
    N_coef[i][0][0] = (x2*y3 - y2*x3)/det; // Must define xi, yi
    N_coef[i][0][1] = (x3*y1 - y3*x1)/det;
    N_coef[i][0][2] = (x1*y2 - y1*x2)/det;
    
    N_coef[i][1][0] = (y2 - y3)/det;
    N_coef[i][1][1] = (y3 - y1)/det;
    N_coef[i][1][2] = (y1 - y2)/det;
    
    N_coef[i][2][0] = (x3 - x2)/det;
    N_coef[i][2][1] = (x1 - x3)/det;
    N_coef[i][2][2] = (x2 - x1)/det;
  }
  
  // Step 4
  for (int i=0; i<M; i++){
    for (int j=0; j<3; j++){
      for (int k=0; k<3; k++){
    
      integral_p = 1; //
      integral_q = 1; //
      integral_r = 1; //
      
      z[i][j][k] = N_coef[i][1][j]*N_coef[i][1][k]*integral_p \
            + N_coef[i][2][j]*N_coef[i][2][k]*integral_q \
            - integral_r;
      }
    
    H[i][j] = -1; //
    }
  }
  
  // Step 5
  for (int i=K+1; i<N; i++){
    for (int j=0; j<3; j++){
      for (int k=0; k<M; k++){
      
        J[i][j][k] = 1; //
      }
    
    I[i][j] = 1; //
    }
  }
  
  // Step 6
  for (int i=0; i<M; i++){
    // Step 7
    for (int k=0; k<3; k++){
      // Step 8
      l = 1; // Find node number for given xk, yk
      
      // Step 9
      if (k > 1){
        for (int j=0; j<k-1; j++){
          // Step 10
          t = 1; // Find node number for given xj, yj
          
          // Step 11
          if (l <= n){
            if (t <= n){
              alpha[l][t] += z[i][k][j];
              alpha[t][l] += z[i][k][j];
            }
            
            else beta[l] -= gamma[t]*z[i][k][j];
          }
          
          else if (t <= n) beta[t] -= gamma[l]*z[i][k][j];
        }
      }
      
      // Step 12
      if (l <= n){
        alpha[l][l] += z[i][k][k];
        beta[l] += H[i][k];
      }
    }
  }
  
  // Step 13
  
  
  
  
}


