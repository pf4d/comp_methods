/*
helper.cpp

Evan M. Cummings
MATH 5310: programming assignment 02
16 March, 2019

Helper functions used for the class.

*/
#include "helper.h"

using namespace std;

// function which outputs the contents of a vector to a file :
void output_vec_to_file(vector<double> & vec,
                        ofstream & file)
{
	for (unsigned int i = 0; i < vec.size(); i++)
	{
		file << std::scientific << vec[i] << '\n';
	}
}

// a function which returns the L^2 norm of a vector :
double norm(vector<double> x)
{
	double val = 0.0;
	for (unsigned int i = 0; i < x.size(); i++)
		val += pow(x[i], 2);
	val = sqrt(val);
	return val;
}



