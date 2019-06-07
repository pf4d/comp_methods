/*
test.cpp

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

Small program which solves the linear system

	\[ A x = b \]

for $x \in \mathhbb{R}^2$ using several iterative methods.

Compile via:

	g++ -o test test.cpp evanmatrix.cpp linearSolver.cpp helper.cpp

*/
#include <vector>
#include "linearSolver.h"
#include "evanmatrix.h"

using namespace std;

// the program :
int main()
{
	/*
	//NOTE: solution is [9.5, 7.5]
	// set the right-hand side vector ``b`` :
	vector<double> b{9, 8};

	// create a matrix with eigenvalues less than 1 :
	vector<double>       vec{0.75, 0.25, 0.25, 0.75};
	vector<unsigned int> rows{0,0,1,1};
	vector<unsigned int> cols{0,1,0,1};
	EvanMatrix A(vec, rows, cols, 2);
	*/

	// three-element poisson problem :
	vector<double>       b{0,1/3.,1/3.,1/6.};
	vector<double>       vec{1,-3,6,-3,-3,6,-3,-3,3};
	vector<unsigned int> rows{0,1,1,1,2,2,2,3,3};
	vector<unsigned int> cols{0,0,1,2,1,2,3,2,3};
	EvanMatrix A(vec, rows, cols, 4);

	A.print();

	// create the solution vector :
	vector<double> x(b.size(), 0.0);

	// show me the diagonal :
	vector<double> diag = A.get_diag();
	for (unsigned int i = 0; i < diag.size(); i++)
		printf(" diag[%u] = %.3e\n", i, diag[i]);

	// type of method to use :
	//string method ("richardson");
	string method ("jacobi");
	//string method ("forwardGaussSeidel");
	//string method ("backwardGaussSeidel");
	//string method ("symmetricGaussSeidel");

	// solve the linear system ``Ax = b`` :
	linear_solve(A, x, b, method, false, 1000, 1e-16, 1e-7);

	// print the solution :
	for (unsigned int i = 0; i < x.size(); i++)
		printf(" x[%i] = %.3e\n", i, x[i]);

	// ensure this is a solution :
	vector<double> b_2 = A.dot(x);

	for (unsigned int i = 0; i < b_2.size(); i++)
		printf(" (b - Ax)_{%u} = %.3e\n", i, b[i] - b_2[i]);

	return 0;
}



