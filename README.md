# Finite Element Method
Parcial 3 de Física Computacional 2 - 2020

## Contributors:
Carolina Herrera Segura, David Silva Sánchez

## Project structure
```
. FiniteElementMethod
├── Makefile
├── examples
|  ├── Integration example
|  └── Linear Equation Solver example
├── include
│  └── headers
├── src
|   └── Finite Element Method Solver
├── tools
|   ├── Integration Methods
|   └── Linear Algebra Methods
├── data
|   ├── input
|   |   ├── data_triangles.txt
|   |   ├── nodes.txt
|   |   └── parameters.txt
|   └── results 
|       ├── gamma_results.txt
|       ├── N_coef_results.txt
|       └── data_eval.txt
└──extras  
    ├── images
    |   └──IntSurface.png
    └── DataVisualization.ipynb
```

## Installation

* Clone the [original repository](https://github.com/DavidSS0397/FiniteElementMethod.git).
* Open a terminal in the cloned repository and run:

```make```

* Run the compiled file ```solver```

```./solver```

## What does it do?

Here we solve a partial differential equation in two dimensions with the general form:

![equation](https://latex.codecogs.com/gif.latex?\frac{\partial}{\partial&space;x}\left(p(x,y)\frac{\partial&space;u}{\partial&space;x}&space;\right)&space;&plus;&space;\frac{\partial}{\partial&space;y}\left(q(x,y)\frac{\partial&space;u}{\partial&space;y}&space;\right)&space;&plus;&space;r(x,y)u(x,y)&space;=&space;f(x,y))

defined in a plane surface D with boundary S. We divide the boundary S into two sub-boundaries S1 and S2. Then, we impose boundary conditions of the form:

On S1:

![equation](https://latex.codecogs.com/gif.latex?u(x,y)=&space;g(x,y))

On S2:

![equation](https://latex.codecogs.com/gif.latex?p(x,y)\frac{\partial&space;u}{\partial&space;x}cos(\theta_1)&plus;q(x,y)\frac{\partial&space;u}{\partial&space;y}cos(\theta_2)&space;&plus;g_1(x,y)u(x,y)&space;=&space;g_2(x,y))

if Φ(x,y) is the approximated solution to u(x,y), the solution given by the algorithm has the form:

![equation](https://latex.codecogs.com/gif.latex?\Phi(x,y)&space;=&space;\sum_{k=0}^{m-1}&space;\gamma_{k}\Phi_k(x,y)&space;\quad&space;(1))

where

![equation](https://latex.codecogs.com/gif.latex?\Phi_k(x,y)&space;=&space;N^{(i)}_j) on T<sub>i</sub> if:

![equation](https://latex.codecogs.com/gif.latex?E_k&space;=&space;(x^{(i)}_j,y^{(i)}_j&space;))

And N<sub>j</sub><sup>(i)</sup> is the polynomial that describes the corresponding edge in the triangle T<sub>i</sub> and has the form:

![equation](https://latex.codecogs.com/gif.latex?N^{(i)}_j(x,y)&space;=&space;a^{(i)}_j&space;&plus;&space;b^{(i)}_jx&space;&plus;&space;c^{(i)}_jy&space;\quad&space;(2))

## How to use it

The code includes a test problem. To use generally, change the files 'data_triangles.txt', 'nodes.txt' and 'parameters.txt' in data/input for the relevant triangle information and change the boundary conditions in src/main.cpp.

The Finite Element Method solves the PDF by dividing the plane of integration in triangular elements and solving the equation for each triangle. The present code doesn't divide the desired plane, so the information of this division must be given as input in the files mentioned above ordered according to the following rules:

* T<sub>0</sub>,...,T<sub>K-1</sub>: internal triangels.
* T<sub>K</sub>,...,T<sub>N-1</sub>: triangles with at least one edge in S2.
* T<sub>N</sub>,...,T<sub>M-1</sub>: triangles with at least one edge in S1.

* E<sub>0</sub>,...,E<sub>p-1</sub>: nodes in S2 that have no boundary with S1.
* E<sub>p</sub>,...,E<sub>n-1</sub>: internal nodes.
* E<sub>n</sub>,...,E<sub>m-1</sub>: nodes in S1.

#### Note
E<sub>i+1</sub> and E<sub>i</sub> have to be consecutive in their corresponding boundary (either S1 or S2). This is important for the line integrals implemented in the method. The order of internal nodes does not matter.

Here is a brief description of the contents of each input data files:

#### nodes.txt

Information relating each node to its x and y coordinates, in that order.

#### data_triangles.txt

Relates the vertices of each triangle to their corresponding nodes and coordinates.

#### parameters.txt

Values for K, N, M, n, p and m, as described above.


Now everything is ready! Compile your code using ```make``` and run ```./solver```.

Your results will be printed in ```data/results/``` .


#### gamma_results.txt

You can see the meaning of gamma in Equation 1. This file contains the values for each entry of the gamma vector.

#### N_coef_results.txt

Coefficients for each polynomial N(x,y) (see Equation 2), ordered by triangle and node/vertex.

#### data_eval.txt

Output of the method ```generate_data```. Matrix of data evaluated using the resulting approximation, given a defined rectangular area.


## Templates for input data

Following this guidelines, here is a brief explanation of how you should fill your nodes.txt and data_triangle.txt files:

#### nodes.txt
```NodeNumber XCOOR YCOOR```

#### data_triangles.txt
```TriangleNumber NodeNumber XCOOR YCOOR```

#### Note
In data_triangles.txt, as each triangle has three vertices, then there are three lines for the same triangle. The order of the nodes is not important, as this is indifferent for the double integrals, just make sure it follows the previously established guidelines.

The NodeNumber, XCOOR and YCOOR in data_triangles.txt and nodes.txt must be equal for each node.

## Example

In order to clarify the last section we will describe the test problem, which is the same as the example proportioned by Burden (section 12.4, see References). The partial differential equation to solve is:

![equation](https://latex.codecogs.com/gif.latex?\frac{\partial^2&space;u(x,y)}{\partial&space;x^2}&space;&plus;&space;\frac{\partial^2&space;u(x,y)}{\partial&space;y^2}&space;=&space;0)

With boundary conditions:

![equation](https://latex.codecogs.com/gif.latex?g(x,y)&space;=&space;4) &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; for (x,y) in L<sub>6</sub> and L<sub>7</sub>

![equation](https://latex.codecogs.com/gif.latex?\frac{\partial&space;u(x,y)}{\partial&space;\mathbf{n}}&space;=&space;x) &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; for (x,y) in L<sub>2</sub> and L4<sub>7</sub>

![equation](https://latex.codecogs.com/gif.latex?\frac{\partial&space;u(x,y)}{\partial&space;\mathbf{n}}&space;=&space;y) &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; for (x,y) in L<sub>6</sub> and L<sub>5</sub>

![equation](https://latex.codecogs.com/gif.latex?\frac{\partial&space;u(x,y)}{\partial&space;\mathbf{n}}&space;=&space;\frac{x&plus;y}{\sqrt2}) &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; for (x,y) in L<sub>1</sub> and L<sub>3</sub>

The following image illustrates how to enumerate the nodes for the given problem. Notice that S2 = L1 + L2 + L4 + L5 and S1 is formed by the rest of the L<sub>i</sub>.

![Integration Surface](https://github.com/DavidSS0397/FiniteElementMethod/blob/master/extras/images/IntSurface.png)

Compare the image with the files in ```data/input``` and the boundary conditions with the functions defined in ```src/main.cpp``` and you'll be ready to test your own system.

### References

The algorithms used in this project where taken from:
Burden, R., & Faires, J. D. (2004). Numerical analysis. Cengage Learning.
