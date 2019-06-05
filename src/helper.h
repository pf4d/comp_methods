/*
helper.cpp

Evan M. Cummings
MATH 5310: programming assignment 02
16 March, 2019

Helper functions used for the class.

*/
#ifndef __HELPER_H
#define __HELPER_H

#include <fstream>
#include <vector>
#include <cmath>

// function which outputs the contents of a vector to a file :
void output_vec_to_file(std::vector<double> & vec, std::ofstream & file);

// a function which returns the L^2 norm of a vector :
double norm(std::vector<double> x);
#endif



