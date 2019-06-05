/*
evanmatrix.h

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

A class representing a sparse matrix in row-major form.

*/
#ifndef __EVANMATRIX_H
#define __EVANMATRIX_H

#include <iostream>
#include <vector>

using namespace std;

class EvanMatrix
{

	public:

		/* constructor */
		EvanMatrix(vector<double>       A,
		           vector<unsigned int> rows,
		           vector<unsigned int> cols,
		           unsigned int         n);

		/* initializer */
		void initialize(vector<double>       A,
		                vector<unsigned int> rows,
		                vector<unsigned int> cols);

		/* matrix dot product function */
		vector<double>       dot(vector<double> & x);

		/* method which prints the non-zero entries */
		void                 print();

		/* accessor functions */
		vector<double>       get_A()     {return A;};
		vector<unsigned int> get_rows()  {return rows;};
		vector<unsigned int> get_cols()  {return cols;};
		unsigned int         get_dim()   {return dim;};
		unsigned int         get_n()     {return n;};
		vector<double>       get_diag()  {return diag;};

	private:

		unsigned int         dim;    // number of non-zero entries
		unsigned int         n;      // dimension of the square matrix
		vector<double>       A;      // non-zero entries of the matrix
		vector<double>       diag;   // n-dimensional diagonal vector
		vector<unsigned int> rows;   // row indices
		vector<unsigned int> cols;   // column indices

		/* method for creating the diagonal vector */
		void                 form_diag();
};
#endif


