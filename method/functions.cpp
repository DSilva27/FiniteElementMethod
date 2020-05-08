#include <iostream>
#include <vector>
#include <fstream>
#include "triangle.h" //??
#include "finite_element.h"
//#include <Armadillo>

using namespace std;
//using namespace arma;

FiniteElement::FiniteElement(){
  cout << "Bienvenido" << endl;
}


FiniteElement::~FiniteElement(){
  cout << "Adios" << endl;
}


void FiniteElement::load_data(){
  
  mat vertices_vec;
  vector <int> nodes_vec;
  vector <int> boundary_vec;
  ifstream file;
  
  int count = 1; // To know in which triangle we are at
  int number, node, bound;
  double x, y, status;
  
  K = 2;
  N = 6;
  M = 10;
  n = 5;
  m = 11;
  
  file.open("../data/data_triangles.txt");
  
  while (file >> number >> node >> x >> y >> bound){
    
    if ( number != count ){
      // When no more data about that triangle, push back element and increase count
      Triangle triangle(vertices_vec, nodes_vec, boundary_vec);
      elements.push_back(triangle);
      
      vertices_vec.clear();
      nodes_vec.clear();
      boundary_vec.clear();
      
      count += 1;
    }
      
    // Fill vectors with data
    vertices_vec.push_back({x, y});
    nodes_vec.push_back(node);
    boundary_vec.push_back(bound);
    
  }
  
  file.close();


  file.open("../data/nodes_v2.txt");

  while (file >> node >> x >> y >> status){

    //Fills the nodes variable with information about each node
    nodes[node].push_back({x, y, status});
  }

  file.close();
}


void FiniteElement::solve(){
  
  // vector <double> gamma(m);
  // vector <double> beta(n);
  // vector <vector <double>> alpha(n, vector <double> (n));
  
  vec gamma(m);
  vec beta(n);
  mat alpha(n, vec(n));
  
  cube N_coef(M, mat(3, vec(3)));
  cube z(M, mat(3, vec(3)));
  cube J((N-K-1), mat(3, vec(3)));
  mat I((N-K-1), vec(3));
  mat H(M, vec (3));
  
  double det;
  double integral_p, integral_q, integral_r;
  
  //double element[i].vertices[0][0], element[i].vertices[2][0], element[i].vertices[2][0];
  //double element[i].vertices[0][1], element[i].vertices[1][1], element[i].vertices[2][1];
  
  int l, t;
  
  
  // Step 1
  for (int l=n+1; l<m; l++){
    gamma[l] = 1; //g(elements[l].vertices[0], vertex[l][0]); // g def is missing

    //    gamma.at(l) = g(vertex[l][0], vertex[l][0]); // g def is missing
  }
  
  // Step 2 is not necessary because vectors are already initialized to 0
  
  // Step 3
  for (int i=0; i<M; i++){
    det = 1; // det function is in LinAlg LinAlg.det(matrix)
    
    // Coefficients of the function N(x, y)
    N_coef[i][0][0] = (elements[i].vertices[1][0]*elements[i].vertices[2][1]\
                       - elements[i].vertices[1][1]*elements[i].vertices[2][0])/det;
    N_coef[i][0][1] = (elements[i].vertices[2][0]*elements[i].vertices[0][1]\
                       - elements[i].vertices[2][1]*elements[i].vertices[0][0])/det;
    N_coef[i][0][2] = (elements[i].vertices[0][0]*elements[i].vertices[1][1]\
                       - elements[i].vertices[0][1]*elements[i].vertices[1][0])/det;
    
    N_coef[i][1][0] = (elements[i].vertices[1][1] - elements[i].vertices[2][1])/det;
    N_coef[i][1][1] = (elements[i].vertices[2][1] - elements[i].vertices[0][1])/det;
    N_coef[i][1][2] = (elements[i].vertices[0][1] - elements[i].vertices[1][1])/det;
    
    N_coef[i][2][0] = (elements[i].vertices[2][0] - elements[i].vertices[1][0])/det;
    N_coef[i][2][1] = (elements[i].vertices[0][0] - elements[i].vertices[2][0])/det;
    N_coef[i][2][2] = (elements[i].vertices[1][0] - elements[i].vertices[0][0])/det;
  }
  
  // Step 4
  for (int i=0; i<M; i++){
    for (int j=0; j<3; j++){
      for (int k=0; k<3; k++){
    
        integral_p = 1; //DblInt( p, elements[i].vertices[0][0], elements[i].vertices[1][0], elements[i].vertices[2][0], 1000, 1000 )
        integral_q = 1; //DblInt( q, elements[i].vertices[0][0], elements[i].vertices[1][0], elements[i].vertices[2][0], 1000, 1000 )
        integral_r = 1; //DblInt( vfunc(r, N_ij, N_ik), elements[i].vertices[0][0], \
                        elements[i].vertices[1][0], elements[i].vertices[2][0], 1000, 1000 )
      
        z[i][j][k] = N_coef[i][1][j]*N_coef[i][1][k]*integral_p \
          + N_coef[i][2][j]*N_coef[i][2][k]*integral_q          \
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
      l = elements[i].nodes[k];
      
      // Step 9
      if (k > 1){
        for (int j=0; j<k-1; j++){
          // Step 10
          t = elements[i].nodes[j];
          
          // Step 11
          if (l <= n){
            if (t <= n){
              alpha[l][t] += z[i][k][j];
              alpha[t][l] += z[i][k][j];

              // alpha.at(l,t) += z[i][k][j];
              // alpha.at(t,l) += z[i][k][j];
            }
            
            else beta[l] -= gamma[t]*z[i][k][j];
            // else beta.at(l) -= gamma.at(t)*z[i][k][j];
          }
          
          else if (t <= n) beta[t] -= gamma[l]*z[i][k][j];
          // else if (t <= n) beta.at(t) -= gamma.at(l)*z[i][k][j];
        }
      }
      
      // Step 12
      if (l <= n){
        alpha[l][l] += z[i][k][k];
        beta[l] += H[i][k];
        
        // alpha.at(l,l) += z[i][k][k];
        // beta.at(l) += H[i][k];
      }
    }
  }
  
  // Step 13
  
  
  
  
}
