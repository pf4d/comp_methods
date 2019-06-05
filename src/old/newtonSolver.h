/*
newtonSolver.cpp

Evan M. Cummings
MATH 5310: programming assignment 02
16 March, 2019

A function which finds roots of a polynomial using Newton-Raphson iteration.

*/
#include <tuple>
#include <vector>
#include <iostream>
#include <cmath>
#include "helper.h"
#include "linearSolver.h"
#include "evanmatrix.h"

using namespace std;

// this is Newton-Raphson polynomial root finding algorithm :
tuple<bool,int> newton_solve(vector<double> f(vector<double> &),
                             EvanMatrix     J(vector<double> &),
                             vector<double> & x,
                             unsigned int   maxiter,
                             double         relaxation_parameter,
                             double         rtol,
                             double         atol)
{
	unsigned int num_dofs          = x.size(); // get num. of dof's
	unsigned int iteration         = 0;        // iteration counter
	bool         converged         = false;    // convergence check
	double       relative_residual = 0.0;      // relative residual of f(x) = 0
	double       residual          = 0.0;      // absolute residual of f(x) = 0
	vector<double> h(num_dofs);                // step size

	// let us know that the algorithm is starting :
	printf("::: beginning non-linear solution process :::\n");

	// begin the solution process :
	while (!converged && iteration < maxiter)
	{
		vector<double> f_i = f(x);  // get the current objective value
		EvanMatrix     J_i = J(x);  // compute the Jacobian

		// perform the linear solve :
		//h = eigen_linear_solve(J_i, f_i);
		h = richardson(J_i, f_i);

		// update solution in the direction `i` :
		for (unsigned int i = 0; i < num_dofs; i++)
		{
			x[i] -= relaxation_parameter * h[i];
		}

		// calculate the absolute residual :
		residual          = norm(f(x));

		// calculate the relative residual :
		relative_residual = norm(h);

		// print statistics to the screen :
   	printf("Newton iteration %d: ",          iteration);
		printf("r (abs) = %.3e (tol = %.3e) ",   residual,          atol);
		printf("r (rel) = %.3e (tol = %.3e)\n",  relative_residual, rtol);

    // increment iteration count :
    iteration++;

		// check for convergence :
		converged = residual < atol || relative_residual < rtol;
	}

	// print the solution when done :
	if (converged)  printf("::: Newton's method converged ");
	else            printf("::: Newton's method did not converge ");

	// always print where the algorithm left off :
	printf("with");
	for (unsigned int i = 0; i < num_dofs; i++)
		printf(" x[%i] = %.3e", i, x[i]);
	printf(" :::\n");

	// finally, return a boolean indicating convergence and the iteration count :
	return make_tuple(converged, iteration);
}



