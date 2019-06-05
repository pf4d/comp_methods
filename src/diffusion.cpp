/*
prob_2.cpp

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

Program which solves the linear system

	\[ A x = b \]

for $x \in \mathhbb{R}^2$ using several iterative methods.

Compile via:

	g++ prob_2.cpp evanmatrix.cpp dirichletbc.cpp -o prob_2.o

*/
#include "diffusion.h"

using namespace std;

double source(double f, double x_l, double x_r)
{
	return f*(x_r - x_l);
}

double stiffness_tensor(unsigned int   i,
                        unsigned int   j,
                        vector<double> k,
                        vector<double> h)
{
	unsigned int n = h.size();
	if (i == 0)
	{
		if      (i == j)           return + k[i] / h[i];
		else if (j == i + 1)       return - k[i] / h[i];
		else                       return 0;
	}
	else if (i > 0 && i < n)
	{
		if      (i == j)           return + k[i-1] / h[i-1] + k[i] / h[i];
		else if (j == i - 1)       return - k[i-1] / h[i-1];
		else if (j == i + 1)       return - k[i] / h[i];
		else                       return 0;
	}
	else
	{
		if      (i == j)           return + k[i-1] / h[i-1];
		else if (j == i - 1)       return - k[i-1] / h[i-1];
		else                       return 0;
	}
}

/* calculates diffusion */
tuple<int,vector<double>> diffusion(Topology       topo,
                                    double         g_r,
                                    double         u_l,
                                    double         f_x,
                                    double         conductivity(double),
                                    string const & solver_method)
{
	double               n   = topo.number_of_vertices();
	double               x_l = topo.x_min();
	vector<double>       x   = topo.coordinates();
	vector<double>       h   = topo.cell_widths();
	vector<double>       u(n);                     // solution vector
	vector<double>       b(n);                     // right-hand side vector
	vector<double>       f(n);                     // source function vector
	vector<double>       k(n);                     // conductivity vector
	vector<double>       K_vals;                   // non-zero values of K tensor
	vector<unsigned int> K_rows;                   // row indices of K_vals
	vector<unsigned int> K_cols;                   // column indices of K_vals
	tuple<bool,int>      linear_out;               // linear solver output


	// next, calculate the cell conductivity :
	for (unsigned int i = 0; i < n; i++)
	{
		double x_mid;
		if      (i == 0)              x_mid = x[i] + h[i] / 2.0;
		else if (i > 0 && i < n-1)    x_mid = x[i] + h[i] / 2.0;
		else                          x_mid = x[i] + h[i] / 2.0;
		k[i] = conductivity(x_mid);
	}

	// next, integrate the source function over elements :
	for (unsigned int i = 0; i < n; i++)
	{
		double x_mid_l, x_mid_r;
		if (i == 0)
		{
			x_mid_l = x[i];
			x_mid_r = x[i] + h[i] / 2.0;
		}
		else
		{
			x_mid_l = x[i-1] + h[i-1] / 2.0;
			x_mid_r = x[i]   + h[i]   / 2.0;
		}
		f[i] = source(f_x, x_mid_l, x_mid_r);
		//printf("i = %u\t x = %f\t k = %f\t f = %f\n", i, x[i], k[i], f[i]);
	}

	// assemble the right-hand side tensor :
	for (unsigned int i = 0; i < n; i++)
	{
		if (i == n-1) b[i] = f[i] - g_r;
		else          b[i] = f[i];
		//printf("i = %u\t b = %f\n", i, b[i]);
	}

	// finally, form the sparse stiffness tensor :
	for (unsigned int i = 0; i < n; i++)
	{
		for (unsigned int j = 0; j < n; j++)
		{
			double K_ij = stiffness_tensor(i, j, k, h);
			if (K_ij != 0)
			{
				K_vals.push_back(K_ij);
				K_rows.push_back(i);
				K_cols.push_back(j);
				//printf("K[%u,%u] = %.3e\n", i, j, K_ij);
			}
		}
	}

	// create a sparse representation of the stiffness tensor :
	EvanMatrix K(K_vals, K_rows, K_cols, n);

	// class which facilitates imposing boundary conditions :
	DirichletBC bc(u_l, x_l, x);

	// apply the boundary condition to the source vector :
	bc.apply(b);

	// apply the boundary condition to the system matrix :
	bc.apply(K);

	/*
	// print the matrix :
	K.print();
	// print the source vector :
	for (unsigned int i = 0; i < b.size(); i++)
		printf(" b[%i] = %.3e\n", i, b[i]);

	// show me the diagonal :
	vector<double> diag = K.get_diag();
	for (unsigned int i = 0; i < diag.size(); i++)
		printf(" diag[%u] = %.3e\n", i, diag[i]);
	*/

	// solve the linear system ``Ax = b`` :
	linear_out = linear_solve(K, u, b,
	                          solver_method,
                            false,
                            1000000000,
                            1e-16,
                            1e-7);

	/*
	// print the solution :
	for (unsigned int i = 0; i < u.size(); i++)
		printf(" u[%i] = %.3e\n", i, u[i]);
	*/

	// calculate the absolute error of the approximation :
	vector<double> b_2 = K.dot(u);
	vector<double> err(u.size());
	for (unsigned int i = 0; i < u.size(); i++)
	{
		err[i] = b[i] - b_2[i];
	}
	double epsilon = norm(err);

	printf("||b - Au|| = %.3e\n", epsilon);

	return make_tuple(get<1>(linear_out), u);
}



