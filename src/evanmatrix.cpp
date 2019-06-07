/*
evanmatrix.cpp

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

A class representing a sparse matrix in row-major form.

*/
#include "evanmatrix.h"

using namespace std;

EvanMatrix::EvanMatrix(vector<double>       A,
                       vector<unsigned int> rows,
                       vector<unsigned int> cols,
                       unsigned int         n)
{
	this->n = n;             // dimension of the square matrix

	// initialize the matrix :
	initialize(A, rows, cols);
}

void EvanMatrix::initialize(vector<double>       A,
                            vector<unsigned int> rows,
                            vector<unsigned int> cols)
{
	this->A    = A;             // non-zero entries of the matrix
	this->rows = rows;          // row indices
	this->cols = cols;          // column indices
	this->dim  = rows.size();   // number of non-zero entries
	this->diag.clear();         // ensure that the diagonal is empty
	this->form_diag();          // create the n-dimensional diagonal vector
}

/* simple function which prints non-zero values of the matrix */
void EvanMatrix::print()
{
	for (unsigned int i = 0; i < A.size(); i++)
	{
		printf("A[%i,%i] = %.2e\n", rows[i], cols[i], A[i]);
	}
}

/* compute the dot product of a vector and this matrix */
vector<double> EvanMatrix::dot(vector<double> & x)
{
	vector<double> b(x.size(), 0.0);  // solution vector

	// compute the next iteration :
	for (unsigned int k = 0; k < this->dim; k++)
	{
		unsigned int i = this->rows[k];    // get the row index
		unsigned int j = this->cols[k];    // get the col index

		b[i] += this->A[k] * x[j];         // that's all there is to it.
	}
	return b;
}

/* create a vector consisting of the diagonal entries of this matrix */
void EvanMatrix::form_diag()
{
	unsigned int i = 0;  // the starting row index

	// iterate through all the non-zero entries of this matrix :
	for (unsigned int k = 0; k < this->dim; k++)
	{
		// if the row index changed without entering a diagonal :
		if (this->rows[k] > i)
		{
			// insert zeros where the diagonal is zero :
			for (unsigned int j = i; j < this->rows[k]; j++)
			{
				diag.push_back(0.0);
				i += 1;
			}
		}
		// if this entry of ``A`` is on the diagonal, insert it to the vector :
		if (this->rows[k] == this->cols[k])
		{
			diag.push_back(this->A[k]);
			i += 1;
		}
	}
	// finally, add on any missing diagonals at the end :
	if (i < this->n)
	{
		// insert zeros where the diagonal is zero :
		for (unsigned int j = i; j < this->n; j++)
		{
			diag.push_back(0.0);
		}
	}
}



