/*
prob_2_part_c.cpp

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

Program which solves the linear system

	\[ A x = b \]

for $x \in \mathhbb{R}^2$ using several iterative methods.

*/
#include <cmath>
#include <vector>
#include <tuple>
#include <ctime>
#include "helper.h"
#include "diffusion.h"
#include "topology.h"

using namespace std;

// conductivity :
double k(double x)
{
	return 1 / (1 - 0.95 * sin(5*M_PI*x));
}

// exact solution :
double u_e(double x)
{
	double val = - 0.5*x*x + x + 19 * cos(5 * M_PI * x) / (100 * M_PI) * (1 - x)
	             + 19 / (500 * M_PI * M_PI) * sin(5 * M_PI * x)
	             + (100 * M_PI - 19) / (100 * M_PI);
	return val;
}

// the program :
int main()
{
	vector<double>            x;              // x-coordinate vector
	vector<double>            u;              // solution vector
	vector<double>            err;            // error vector
	tuple<int,vector<double>> diff_out;       // diffusion output tuple

	double epsilon;                           // norm of error
	unsigned int n;                           // number of dofs

	clock_t t_0, t_f;                         // for the timer

	// the num. of divisions of ``x`` to create
	vector<double> n_a{10, 20, 40, 80, 160, 320, 640};

	// initialize the error and iteration vector per refinement :
	vector<double> err_a(n_a.size());
	vector<double> itr_a(n_a.size());
	vector<double> tme_a(n_a.size());

	// type of method to use :
	string solver_method = "forwardGaussSeidel";

	// iterate over each refinement :
	for (unsigned int i = 0; i < n_a.size(); i++)
	{
		// get the number of divisions :
		n = n_a[i];

		double       x_l = 0.0;         // left endpoint
		double       x_r = 1.0;         // right endpoint
		double       g_r = 0.0;         // outward flux of u at x = x_r
		double       u_l = 1.0;         // left essential boundary condition
		double       f_x = 1.0;         // constant source function

		// create the topology for the domain described above :
		Topology topo(n, x_l, x_r);

		// retrieve the vector of x coordinates :
		vector<double> x = topo.coordinates();

		// solve the diffusion problem :
		t_0 = clock();
		diff_out = diffusion(topo, g_r, u_l, f_x, k, solver_method);
		t_f = clock();

		// extract the solution vector :
		u = get<1>(diff_out);

		// calculate the error of the approximation :
		err.resize(n);
		for (unsigned int j = 0; j < n; j++)
		{
			err[j] = u_e(x[j]) - u[j];
		}

		// calculate the norm of the error :
		epsilon = norm(err) / n;

		// print the error :
		printf("||u_e - u|| = %.3e\n", epsilon);

		// append the error to the vector :
		err_a[i] = epsilon;
		itr_a[i] = get<0>(diff_out);
		tme_a[i] = (double) (t_f - t_0) / CLOCKS_PER_SEC;

		// output the vectors to files for later plotting :
		ofstream x_file ("../data/c/x_" + to_string(n) + ".txt");
		ofstream u_file ("../data/c/u_" + to_string(n) + ".txt");

		output_vec_to_file(x, x_file);
		output_vec_to_file(u, u_file);

		x_file.close();
		u_file.close();
	}
	// output the vectors to files for later plotting :
	ofstream n_file ("../data/c/n.txt");
	ofstream e_file ("../data/c/e.txt");
	ofstream i_file ("../data/c/i.txt");
	ofstream t_file ("../data/c/t.txt");

	output_vec_to_file(n_a,   n_file);
	output_vec_to_file(err_a, e_file);
	output_vec_to_file(itr_a, i_file);
	output_vec_to_file(tme_a, t_file);

	n_file.close();
	e_file.close();
	i_file.close();
	t_file.close();

	return 0;
}



