/*
evanmatrix.h

Evan M. Cummings
MATH 5310: programming assignment 03
9 May, 2019

A class representing a topological space.

*/
#include <vector>
#include "topology.h"

using namespace std;

/* constructor */
Topology::Topology(unsigned int n,
                   double       x_l,
                   double       x_r)
{
	this->n   = n;           // number of vertices
	this->x_l = x_l;         // left endpoint
	this->x_r = x_r;         // right endpoint

	this->n_e = n - 1;       // number of elements
	this->l   = x_r - x_l;   // width of the domain
	this->h_x = l / n_e;     // width of individual cells (constant for now)

	this->x.resize(n);       // resize the vertex coordinate vector
	this->h.resize(n_e);     // resize the cell width vector

	// first, initialize the cell widths (constant for now) :
	for (unsigned int i = 0; i < n_e; i++)
	{
		h[i] = h_x;
	}

	// then the x-coordinate of each degree of freedom :
	for (unsigned int i = 0; i < n; i++)
	{
		if (i == 0) x[i] = x_l;
		else        x[i] = x[i-1] + h[i-1];
	}

};



