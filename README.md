# Finite Element Method
Parcial 3 de Física Computacional 2 - 2020

### Contributors:
Carolina Herrera Segura, David Silva Sánchez

### Structure of the project
```
. FiniteElementMethod
├── Makefile
├── examples
|  ├── Integration example
|  └── Linear Equation Solver example
├── include
│  └── headers
└── src
|   └── Finite Element Method Solver
└── tools
|   ├── Integration Methods
|   └── Linear Algebra Methods
└── data
    ├── Input
    |   ├── data_triangles.txt
    |   └── nodes.txt
    └── Output    
    |   ├── gamma_results.txt
    |   └── N_coef_results.txt
```

### Instalation

* Clone the [original repository](https://github.com/DavidSS0397/FiniteElementMethod.git).
* Open a terminal where you cloned the repository.
* Now follow the following guidelines on how to use it.

### What does it do?

Here we solve a partial differential equation in two dimensions with the general expression:

![equation](https://latex.codecogs.com/gif.latex?\frac{\partial}{\partial&space;x}\left(p(x,y)\frac{\partial&space;u}{\partial&space;x}&space;\right)&space;&plus;&space;\frac{\partial}{\partial&space;y}\left(q(x,y)\frac{\partial&space;u}{\partial&space;y}&space;\right)&space;&plus;&space;r(x,y)u(x,y)&space;=&space;f(x,y))

If S is the boundary of a plane surface D, we divide the boundary S into two sub-boundaries S1 and S2. Then, the function u(x,y) must have boundary conditions of the form:

In S1:

![equation](https://latex.codecogs.com/gif.latex?u(x,y)=&space;g(x,y))

In S2:

![equation](https://latex.codecogs.com/gif.latex?p(x,y)\frac{\partial&space;u}{\partial&space;x}cos(\theta_1)&plus;q(x,y)\frac{\partial&space;u}{\partial&space;y}cos(\theta_2)&space;&plus;g_1(x,y)u(x,y)&space;=&space;g_2(x,y))

if Φ(x,y) is the approximated solution to u(x,y), the solution given by the algorithm has the form:

![equation](https://latex.codecogs.com/gif.latex?\Phi(x,y)&space;=&space;\sum_{k=0}^{m-1}&space;\gamma_{k}\Phi_k(x,y)&space;\quad&space;(1))

where

![equation](https://latex.codecogs.com/gif.latex?\Phi_k(x,y)&space;=&space;N^{(i)}_j) on T<sub>i</sub> if:

![equation](https://latex.codecogs.com/gif.latex?E_k&space;=&space;(x^{(i)}_j,y^{(i)}_j&space;))

And N<sup>(i)</sup><sub>j</sub> is the polynomial that describes the corresponding edge in the triangle T<sub>i</sub> and has the form:

![equation](https://latex.codecogs.com/gif.latex?N^{(i)}_j(x,y)&space;=&space;a^{(i)}_j&space;&plus;&space;b^{(i)}_jx&space;&plus;&space;c^{(i)}_jy&space;\quad&space;(2))

When you use the code, you have to include the functions for your problem. This has to be done in src/main.cpp. Here you just have to change the definitions of each function.

The Finite Element Method solves the PDF by dividing the plane of integration in triangles and solving the equation for each triangles. This algorithm doesn't divide the desired plane, this has to be done manually and the information has to be ordered according to the following rules. Templates are found in data/input .

* T<sub>0</sub>,...,T<sub>K-1</sub>: number of internal triangels
* T<sub>K</sub>,...,T<sub>N-1</sub>: number of triangles with at least one edge in S2
* T<sub>N</sub>,...,T<sub>M-1</sub>: number of triangles with at least one edge in S1

* E<sub>0</sub>,...,E<sub>p-1</sub>: number of nodes in S2 that have no boundary with S1
* E<sub>p</sub>,...,E<sub>n-1</sub>: number of nodes internal nodes
* E<sub>n</sub>,...,E<sub>m-1</sub>: number of nodes in S1

Note: E<sub>i+1</sub> and E<sub>i</sub> have to be consecutive in their corresponding boundary (either S1 or S2). This is important for the line integrals implemented in the method (this doesn't apply for internal nodes).

Following this guidelines, here is a brief explanation of how you should fill your nodes.txt and data_triangle.txt files:

#### nodes.txt
```NodeNumber XNode YNode}```

#### data_triangles.txt
```TriangleNumber NodeNumber XNode YNode```

Note: in data_triangles.txt as each triangle has three vertices, then there are three lines for the same triangles. The order in which you input the nodes in this file is not important, as this is indifferent for the double integrales. Just make sure it follows the guidelines stablished previously. 

The NodeNumber, XNode and YNode in data_triangles.txt and nodes.txt have to coincide.

Now everything is ready. Compile your code using ```make``` and run ```./solver```.

Your results will be printed in data/results/ .

### gamma_results.txt

You can see the meaning of γ in Equation 1. This file has the values of each value of the gamma vector.

### N_coef_results.txt

This file has the coefficients for each polynome N (see Equation 2).


### Acknowledgments

The algorithms used in this project where taken from:
Burden, R., & Faires, J. D. (2004). Numerical analysis. Cengage Learning.
