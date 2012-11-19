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

#ifndef NEWTONS_METHOD_HPP
#define NEWTONS_METHOD_HPP

#include <vector>
#include <numeric>


// unrestricted newtons method
bool
newtons_method(const bool VERBOSE,
	       const std::vector<double> &initial_lambdas,
	       const std::vector<double> &initial_xs,
	       const std::vector<double> &moments,
	       const double values_sum,
	       const double tolerance, const size_t max_iter,
	       std::vector<double> &root_lambdas,
	       std::vector<double> &root_xs);


// newton's method so that all elements are strictly positive

#endif
