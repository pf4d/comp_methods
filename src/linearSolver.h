/*
linearSolver.cpp

Evan M. Cummings
MATH 5310: programming assignment 03
1 May, 2019

A function which solves the linear system ``A x = b`` for ``x``.

*/
#ifndef __LINEARSOLVER_H
#define __LINEARSOLVER_H

#include <vector>
#include <tuple>
#include <cmath>
#include "evanmatrix.h"
#include "helper.h"

using namespace std;

// general iterative method for solving linear systems of equations :
tuple<bool,int> linear_solve(EvanMatrix     & A,
                             vector<double> & x,
                             vector<double> & b,
                             string const   & method,
                             bool             verbose    = false,
                             unsigned int     max_iter   = 1000,
                             double           atol       = 1e-16,
                             double           rtol       = 1e-7);
#endif



