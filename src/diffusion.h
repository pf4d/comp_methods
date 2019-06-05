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
#ifndef __DIFFUSION_H
#define __DIFFUSION_H

#include <cmath>
#include <vector>
#include "helper.h"
#include "linearSolver.h"
#include "dirichletbc.h"
#include "evanmatrix.h"
#include "topology.h"

using namespace std;

double source(double f, double x_l, double x_r);

double stiffness_tensor(unsigned int   i,
                        unsigned int   j,
                        vector<double> k,
                        vector<double> h);

/* calculates diffusion */
tuple<int,vector<double>> diffusion(Topology       topo,
                                    double         g_r,
                                    double         u_l,
                                    double         f_x,
                                    double         conductivity(double),
                                    string const & solver_method);
#endif



