#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

//#include <dolfin/common/Array.h>

#include "helper.h"
#include "evanmatrix.h"
#include "linearSolver.h"
#include "diffusion.h"

namespace py = pybind11;

PYBIND11_MODULE(comp_methods, m)
{
	m.def("linear_solve", & linear_solve);

	m.def("diffusion", & diffusion);

	m.def("output_vec_to_file", & output_vec_to_file);

	//FIXME: this doesn't work :
	//m.def("norm", & norm);

	py::class_<EvanMatrix>(m, "EvanMatrix")
		.def(py::init<vector<double>,
		           vector<unsigned int>,
		           vector<unsigned int>,
		           unsigned int>())
		.def("dot",        &EvanMatrix::dot)
		.def("get_diag",   &EvanMatrix::get_diag)
		.def("print",      &EvanMatrix::print)
		.def("initialize", &EvanMatrix::initialize)
		.def("get_A",      &EvanMatrix::get_A)
		.def("get_rows",   &EvanMatrix::get_rows)
		.def("get_cols",   &EvanMatrix::get_cols)
		.def("get_dim",    &EvanMatrix::get_dim)
		.def("get_n",      &EvanMatrix::get_n);

	py::class_<DirichletBC>(m, "DirichletBC")
		.def(py::init<double, double, vector<double> &>())
		.def("get_val",      &DirichletBC::get_val)
		.def("eval",         &DirichletBC::eval)
		.def("form_indices", &DirichletBC::form_indices)
		.def("apply",
		     (void (DirichletBC::*)(vector<double> &)) &DirichletBC::apply)
		.def("apply",
		     (void (DirichletBC::*)(EvanMatrix &)) &DirichletBC::apply);

	py::class_<Topology>(m, "Topology")
		.def(py::init<unsigned int, double, double>())
		.def("number_of_elements", &Topology::number_of_elements)
		.def("number_of_vertices", &Topology::number_of_vertices)
		.def("coordinates",        &Topology::coordinates)
		.def("cell_widths",        &Topology::cell_widths)
		.def("x_min",              &Topology::x_min)
		.def("x_max",              &Topology::x_max)
		.def("width",              &Topology::width);
}



