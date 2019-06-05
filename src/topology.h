/*
evanmatrix.h

Evan M. Cummings
MATH 5310: programming assignment 03
9 May, 2019

A class representing a topological space.

*/
#ifndef __TOPOLOGY_H
#define __TOPOLOGY_H

#include <vector>

using namespace std;

class Topology
{

	public:

		/* constructor */
		Topology(unsigned int n,
		         double       x_l,
		         double       x_r);

		/* accessor functions */
		unsigned int   number_of_elements()  {return n_e;};
		unsigned int   number_of_vertices()  {return n;};
		vector<double> coordinates()         {return x;};
		vector<double> cell_widths()         {return h;};
		double         x_min()               {return x_l;};
		double         x_max()               {return x_r;};
		double         width()               {return l;};

	private:

		unsigned int    n;         // number of vertices
		unsigned int    n_e;       // number of elements
		double          h_x;       // width of individual cells
		double          x_r;       // right endpoint
		double          x_l;       // left endpoint
		double          l;         // domain width

		vector<double>  x;         // x-coordinate vector
		vector<double>  h;         // cell-width vector

};
#endif



