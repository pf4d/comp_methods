/*
dirichletbc.h

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

A class representing an essential boundary condition.

*/
#ifndef __DIRICHLETBC_H
#define __DIRICHLETBC_H

#include <vector>
#include <cmath>
#include "evanmatrix.h"

using namespace std;

class DirichletBC
{

	public:

		/* constructor */
		DirichletBC(double           val,
                double           x_0,
                vector<double> & x);

		/* accessor function */
		double get_val()   {return val;};

		/* evaluation function */
		bool   eval(double x, double tol = 1e-10);

		/* function which collects indicies at boundaries */
		void   form_indices();

		/* application function */
		void   apply(vector<double> & u);
		void   apply(EvanMatrix     & A);

	private:

		double               val;
		double               x_0;
		vector<double>       x;
		vector<unsigned int> idx;

};
#endif



