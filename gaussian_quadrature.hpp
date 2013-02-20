/*    Copyright (C) 2012 University of Southern California and
 *                       Andrew D. Smith and Timothy Daley
 *
 *    Authors: Andrew D. Smith and Timothy Daley
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GAUSSIAN_QUADRATURE_HPP
#define GAUSSIAN_QUADRATURE_HPP

#include <vector>
#include <numeric>

// Andrew's code for gaussian quadrature
// finds the points of weight by finding
// roots of the associated polynomial
void poly_solve_gauss_quad(const std::vector<double> &moments,
			   const size_t n_points,
			   std::vector<double> &weights,
			   std::vector<double> &points);

void
golub_welsh_quadrature(const std::vector<double> &moments,
		       const size_t n_points,
		       const double tol, const size_t max_iter,
		       std::vector<double> &points,
		       std::vector<double> &weights);

void
laguerre_modified_quadrature(const std::vector<double> &moments,
			     const size_t n_points,
			     const double tol, const size_t max_iter,
			     std::vector<double> &points,
			     std::vector<double> &weights);



#endif
