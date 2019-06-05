/*
newtonSolver.cpp

Evan M. Cummings
MATH 5310: programming assignment 02
16 March, 2019

A function which finds roots of a polynomial using Newton-Raphson iteration.

*/
#ifndef __NEWTONSOLVER_H
#define __NEWTONSOLVER_H

using namespace std;

// this is Newton-Raphson polynomial root finding algorithm :
tuple<bool,int> newton_solve(vector<double> f(vector<double> &),
                             EvanMatrix     J(vector<double> &),
                             vector<double> & x,
                             unsigned int   maxiter,
                             double         relaxation_parameter,
                             double         rtol,
                             double         atol);
#endif


