/*
prob_1.cpp

Evan M. Cummings
MATH 5310: programming assignment 02
21 April, 2019

Small program which finds roots of the polynomial

	\[ f(z) = 1 + z^2 \exp(z) \]

for $z \in \mathhbb{C}$ using Newton-Raphson iteration.

Compile via:

	g++ -I /path/to/eigen3 prob_1.cpp -o prob_1.o
	g++ -I /home/pf4d/local/eigen-3.3.7/include/eigen3/ prob_1.cpp evanmatrix.cpp -o prob_1.o

*/
#include <cmath>
#include <vector>
#include <tuple>
#include "helper.h"
#include "newtonSolver.h"
#include "evanmatrix.h"

using namespace std;

// the scalar function ``F_1`` :
double F_1(vector<double> & x)
{
	double val = 1 + pow(x[0], 2) - pow(x[1], 2) + exp(x[0])*cos(x[1]);
	return val;
}

// the scalar function ``F_2`` :
double F_2(vector<double> & x)
{
	double val = 2*x[0]*x[1] + exp(x[0])*sin(x[1]);
	return val;
}

// the (scalar) partial derivative of ``F_1`` with respect to ``x`` :
double dF_1_dx(vector<double> & x)
{
	double val = 2*x[0] + exp(x[0])*cos(x[1]);
	return val;
}

// the (scalar) partial derivative of ``F_2`` with respect to ``x`` :
double dF_2_dx(vector<double> & x)
{
	double val = 2*x[1] + exp(x[0])*sin(x[1]);
	return val;
}

// the vector function `F` :
vector<double> F(vector<double> & x)
{
	vector<double> vec{F_1(x), F_2(x)};
	return vec;
}

// the Jacobian of `F` :
// NOTE: return matrix as a vector column-major format
EvanMatrix J(vector<double> & x)
{
	double J_11 = dF_1_dx(x);
	double J_21 = dF_2_dx(x);
	vector<double> vec{J_11, -J_21, J_21, J_11};
	vector<unsigned int> rows{0,0,1,1};
	vector<unsigned int> cols{0,1,0,1};
	EvanMatrix mat(vec, rows, cols);
	return mat;
}

// the program :
int main()
{
	// set the initial condition ``z`` :
	vector<double> z{-1, 4};

	// solve for the zero of `f(z)` :
	newton_solve(F, J, z, 2, 1.0, 1e-10, 1e-16);

	return 0;
}



