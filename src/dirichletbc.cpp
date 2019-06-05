/*
dirichletbc.cpp

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

A class representing an essential boundary condition.

*/

#include "dirichletbc.h"

using namespace std;

/* constructor */
DirichletBC::DirichletBC(double           val,
                         double           x_0,
                         vector<double> & x)
{
	this->val = val;       // the value at the boundary
	this->x_0 = x_0;       // the coordinate of the boundary
	this->x   = x;         // the vector of coordinates

	this->form_indices();  // where all the magic happens
}

/* this method really isn't needed, but might later */
bool DirichletBC::eval(double x, double tol)
{
	if (fabs(x - this->x_0) < tol)  return true;
	else                            return false;
}

/* method applies the boundary condition to a vector */
void DirichletBC::apply(vector<double> & u)
{
	// iterate through each index and set the vector equal to the desired value :
	for (unsigned int i = 0; i < this->idx.size(); i++)
	{
		u[this->idx[i]] = this->val;
	}
}

/* method applies the boundary condition to a tensor */
void DirichletBC::apply(EvanMatrix & A)
{
	unsigned int         dim     = A.get_dim();    // number of non zeros
	vector<double>       vals    = A.get_A();      // non-zero values
	vector<unsigned int> rows    = A.get_rows();   // corresponding row indices
	vector<unsigned int> cols    = A.get_cols();   // corresponding col indices

	// iterate through all the indicies that have been marked :
	for (unsigned int i = 0; i < this->idx.size(); i++)
	{
		unsigned int k = this->idx[i];  // nodal boundary markers

		// iterate though all the non-zero entries of the tensor :
		for (unsigned int j = 0; j < dim; j++)
		{
			// if the entry is on the diagonal, make it equal ``1`` :
			if      (rows[j] == k && cols[j] == k)
			{
				vals[j] = 1.0;
			}
			// otherwise it is zero, so erase it from the sparse matrix :
			else if (rows[j] == k && cols[j] != k)
			{
				vals.erase(vals.begin() + j);
				rows.erase(rows.begin() + j);
				cols.erase(cols.begin() + j);
			}
		}
	}
	// finally, re-initialize the matrix with new values :
	A.initialize(vals, rows, cols);
}

/* form a vector of dirichlet nodes */
void DirichletBC::form_indices()
{
	// iterate over each dof :
	for (unsigned int i = 0; i < this->x.size(); i++)
	{
		// is this on the boundary ?
		if (eval(this->x[i]))
		{
			this->idx.push_back(i);
		}
	}
}



