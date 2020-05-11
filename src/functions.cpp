#include <iostream>
#include <vector>
#include <fstream>
#include "triangle.h"
#include "finite_element.h"
#include "DblInt.h"
#include "LInt.h"
#include "LinAlg.h"

using namespace std;

FiniteElement::FiniteElement(){
  //cout << "Bienvenido" << endl;
}


FiniteElement::~FiniteElement(){
  //cout << "Adios" << endl;
}


void FiniteElement::load_data(){
  
  mat vertices_vec;
  vector <int> nodes_vec;
  ifstream file;
  
  int number, node;
  double x, y;
  
  K = 2;
  N = 6;
  M = 10;
  n = 5;
  m = 11;
  p = 4;
  
  // Save data about each triangle to Triangle objects
  file.open("data/input/data_triangles.txt");
  
  for (int i=0; i<M; i++){//
  
    for (int j=0; j<3; j++){
      // Fill vectors with data
      file >> number >> node >> x >> y;
      
      vertices_vec.push_back({x, y});
      nodes_vec.push_back(node);
    }
    
    // Create Triangle object with data
    Triangle triangle(vertices_vec, nodes_vec);
    elements.push_back(triangle);
    
    vertices_vec.clear();
    nodes_vec.clear();
  }
  
  file.close();
  
  // Save coordinates of each node to matrix
  file.open("data/input/nodes.txt");
  
  while (file >> node >> x >> y){
  
    //Fills nodes vector with information about each node
    nodes.push_back({x, y});
  }
  
  file.close();
}


void FiniteElement::solve( vfunc VF){
  
  vec gamma(m);
  vec beta(n);
  mat alpha(n, vec(n));
  
  cube N_coef(M, mat(3, vec(3)));
  cube z(M, mat(3, vec(3)));
  cube J((N-K), mat(3, vec(3)));
  mat I((N-K), vec(3));
  mat H(M, vec (3));
  
  double det;
  double integral_p, integral_q, integral_r;
  
  int step = 200;
  int l, t;
  
  ofstream gamma_file( "data/results/gamma_results.txt", ios::out );
  ofstream N_coef_file( "data/results/N_coef_results.txt", ios::out );
  
  if ( !gamma_file or !N_coef_file) 
    {
      cout << "Can't open file'" << endl;
      exit( 1 );
    }
  
  // Step 1
  for (int l=n; l<m; l++){
    gamma[l] = VF[4](nodes[l][0], nodes[l][1]);
  }
  
  // Step 2 is not necessary because vectors are already initialized to 0
  
  // Step 3
  for (int i=0; i<M; i++){

    mat Matrix{ {1, elements[i].vertices[0][0], elements[i].vertices[0][1]},
                {1, elements[i].vertices[1][0], elements[i].vertices[1][1]},
                {1, elements[i].vertices[2][0], elements[i].vertices[2][1]} };

    det = LinAlg.Det33(Matrix); 
    
    // Coefficients of the function N(x, y)
    // N[i][a,b,c][1,2,3]
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
      for (int k=0; k<=j; k++){
        
        integral_p = DInt.DoubleIntegral( VF[0], elements[i].vertices, step, step );
        integral_q = DInt.DoubleIntegral( VF[1], elements[i].vertices, step, step );
        integral_r = DInt.DoubleIntegral( VF[2], N_coef[i][j], N_coef[i][k], elements[i].vertices, step, step );
      
        z[i][j][k] = N_coef[i][1][j]*N_coef[i][1][k]*integral_p \
                    + N_coef[i][2][j]*N_coef[i][2][k]*integral_q \
                    - integral_r;
      }
      
      H[i][j] = -DInt.DoubleIntegral( VF[3], N_coef[i][j], elements[i].vertices, step, step ); //
    }
  }
  
  // Step 5
  for (int i=K; i<N; i++){
    for (int j=0; j<3; j++){
      for (int k=0; k<=j; k++){
        
        double Int = 0;
        
        for (int l=0; l<p; l++){
          
          if (l == 0){
            Int += LInt.LineIntegral(VF[5], N_coef[i][j], N_coef[i][k], nodes[n], nodes[l], step);
          }
          
          else if(l == p-1){
            Int += LInt.LineIntegral(VF[5], N_coef[i][j], N_coef[i][k], nodes[l], nodes[m-1], step);
          }
          
          else{
            Int += LInt.LineIntegral(VF[5], N_coef[i][j], N_coef[i][k], nodes[l], nodes[l+1], step);
          }
        }
        
        J[i-K][j][k] = Int;
      }
     
      double Int = 0;
      
      for (int l = 0; l < n; l++){
        if (l == 0){
          Int += LInt.LineIntegral( VF[6], N_coef[i][j], nodes[n], nodes[l], step);
        }
        
        else if(l == n-1){
          Int += LInt.LineIntegral( VF[6], N_coef[i][j], nodes[l], nodes[m-1], step);
        }
        
        else{
          Int += LInt.LineIntegral( VF[6], N_coef[i][j], nodes[l], nodes[l+1], step);
        }
      }
      
      I[i-K][j] = Int;
    }
  }
  
  // Step 6
  for (int i=0; i<M; i++){
    // Step 7
    for (int k=0; k<3; k++){
      // Step 8
      l = elements[i].nodes[k];
      
      // Step 9
      if (k > 0){
        for (int j=0; j<k; j++){
          // Step 10
          t = elements[i].nodes[j];
          
          // Step 11
          if (l < n){
            if (t < n){
              alpha[l][t] += z[i][k][j];
              alpha[t][l] += z[i][k][j];
            }
            
            else beta[l] -= gamma[t]*z[i][k][j];
          }
          
          else if (t < n) beta[t] -= gamma[l]*z[i][k][j];
        }
      }
      
      // Step 12
      if (l < n){
        alpha[l][l] += z[i][k][k];
        beta[l] += H[i][k];
      }
    }
  }
  
  // Step 13
  for (int i=K; i<N; i++ ){
  // Step 14
    for (int k=0; k<3; k++){
      //Step 15
      l = elements[i].nodes[k];
      
      //Step 16
      if (k > 0){
        for (int j=0; j<k-1; j++){
          // Step 17
          t = elements[i].nodes[j];
          
          // Step 18
          if (l < n){
            if (t < n){
              alpha[l][t] += J[i-K][k][j];
              alpha[t][l] += J[i-K][k][j];
            }
            
            else beta[l] -= gamma[t]*J[i-K][k][j];
          }
          
          else if (t < n) beta[t] -= gamma[l]*J[i-K][k][j];
        }
      }
      
      // Step 19
        if (l < n){
          alpha[l][l] += J[i-K][k][k];
          beta[l] += I[i-K][k];
        }
    }
  }
  // Step 20
  LinAlg.SOR(alpha, beta, gamma, 1.25, 0.00003, 2000);
  
  // Step 21
  // Save gamma entries to file
  for (int i=0; i<m; i++){
    gamma_file << gamma[i] << " ";
    }
  gamma_file << endl;
  
  gamma_file.close();
  
  // Save N coefficients to file (number, ai, bi, ci)
  N_coef_file << "TRIANGLE ai bi ci" << endl;
  
  for (int i=0; i<M; i++){
    for (int j=0; j<3; j++){
      N_coef_file << i;
      for (int k=0; k<3; k++){
        N_coef_file << " " << N_coef[i][k][j];
      }
      
      N_coef_file << endl;
    }
  }
  
  N_coef_file.close();
}
