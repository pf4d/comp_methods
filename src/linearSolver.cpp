/*
linearSolver.cpp

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

A function which solves the linear system ``A x = b`` for ``x``.

*/
#include "linearSolver.h"

using namespace std;

// general iterative method for solving linear systems of equations :
tuple<bool,int> linear_solve(EvanMatrix     & A,
                             vector<double> & x,
                             vector<double> & b,
                             string const   & method,
                             bool             verbose,
                             unsigned int     max_iter,
                             double           atol,
                             double           rtol)
{
	unsigned int dim            = b.size();     // number of dof
	unsigned int num_non_zero   = A.get_dim();  // the number of non-zeros of A
	unsigned int iteration      = 0;            // iteration count
	bool         converged      = false;        // convergence flag

	double       residual;                      // iteration residual
	double       residual_0;                    // initial residual
	double       relative_residual;             // relative residual

	vector<double>       A_vals = A.get_A();    // matrix value vector
	vector<double>       A_diag = A.get_diag(); // the diagonal entries of A
	vector<unsigned int> A_rows = A.get_rows(); // row indices of A_vals
	vector<unsigned int> A_cols = A.get_cols(); // col indices of A_vals

	vector<double> x_n(dim, 0.0);               // new solution vector
	vector<double> x_s(dim, 0.0);               // intermediate solution vector
	//vector<double> x(dim,   0.0);               // prev. sol'n vector
	vector<double> r(dim,   0.0);               // residual vector.
	vector<double> A_dot_x_n(dim);              // dot product of A with x_n

	fill(x.begin(), x.end(), 0.0);              // initialize the guess to zero

	unsigned int i,j,k,k_b,i_b,j_b;             // various loop indices

	// let us know that the algorithm is starting :
	printf("::: solving %u x %u linear system using %s iteration :::\n",
	       dim, dim, method.c_str());

	// begin the solution process :
	while (!converged && iteration < max_iter)
	{

		// reset the current solution vectors :
		for (i = 0; i < dim; i++)
		{
			x_n[i] = 0.0;
			x_s[i] = 0.0;
		}

		// compute the next iteration :
		for (k = 0; k < num_non_zero; k++)
		{
			i = A_rows[k];    // get the current row index
			j = A_cols[k];    // get the current col index

			// calc. reversed indices if using backward or symmetric Gauss-Seidel :
			if (method == "backwardGaussSeidel" || method == "symmetricGaussSeidel")
			{
				k_b = num_non_zero - 1 - k;
				i_b = A_rows[k_b];    // get the current row index
				j_b = A_cols[k_b];    // get the current col index
			}

			// calculate the next iterate :
			if      (method == "richardson")
			{
				if (i == j)     x_n[i] += (1 - A_vals[k])*x[i] + b[i];
				else            x_n[i] -= A_vals[k] * x[j];
			}
			else if (method == "jacobi")
			{
				if (i == j)     x_n[i] += b[i] / A_diag[i];
				else            x_n[i] -= A_vals[k] / A_diag[i] * x[j];
			}
			else if (method == "forwardGaussSeidel")
			{
				if (i == j)     x_n[i] += b[i] / A_diag[i];
				else if (j < i) x_n[i] -= A_vals[k] / A_diag[i] * x_n[j];
				else if (j > i) x_n[i] -= A_vals[k] / A_diag[i] * x[j];
			}
			else if (method == "backwardGaussSeidel")
			{
				if (i_b == j_b)     x_n[i_b] += b[i_b] / A_diag[i_b];
				else if (j_b < i_b) x_n[i_b] -= A_vals[k_b] / A_diag[i_b] * x[j_b];
				else if (j_b > i_b) x_n[i_b] -= A_vals[k_b] / A_diag[i_b] * x_n[j_b];
			}
			else if (method == "symmetricGaussSeidel")
			{
				if (i == j)     x_s[i] += b[i] / A_diag[i];
				else if (j < i) x_s[i] -= A_vals[k] / A_diag[i] * x_n[j];
				else if (j > i) x_s[i] -= A_vals[k] / A_diag[i] * x[j];

				if (i_b == j_b)     x_n[i_b] += b[i_b] / A_diag[i_b];
				else if (j_b < i_b) x_n[i_b] -= A_vals[k_b] / A_diag[i_b] * x_s[j_b];
				else if (j_b > i_b) x_n[i_b] -= A_vals[k_b] / A_diag[i_b] * x_n[j_b];
			}
		}

		// calculate the absolute residual vector :
		A_dot_x_n = A.dot(x_n);
		for (i = 0; i < dim; i++)
		{
			r[i] = b[i] - A_dot_x_n[i];
		}

		// calculate the absolute residual :
		residual = norm(r);

		// save the first residual :
		if (iteration == 0)    residual_0 = residual;

		// calculate the relative residual :
		relative_residual = residual / residual_0;

		// print statistics to the screen :
		if (verbose)
		{
			printf("iteration %d: ",                 iteration);
			printf("r (abs) = %.3e (tol = %.3e) ",   residual,          atol);
			printf("r (rel) = %.3e (tol = %.3e)\n",  relative_residual, rtol);
		}

		// set the previous solution for the next iteration :
		for (i = 0; i < dim; i++)
		{
			x[i] = x_n[i];
		}

    // increment iteration count :
    iteration++;

		// check for convergence :
		converged = residual < atol || relative_residual < rtol;
	}

	// print the solution when done :
	if (converged)
		printf("::: Iterative method converged in %u iterations :::\n",
		       iteration-1);
	else
		printf("::: Iterative method did not converge :::\n");

	// return the solution :
	return make_tuple(converged, iteration);
}



