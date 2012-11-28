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

#ifndef LIBRARY_SIZE_ESTIMATES_HPP
#define LIBRARY_SIZE_ESTIMATES_HPP

#include <vector>
#include <numeric>

#include "newtons_method.hpp"

double
chao87_lowerbound_librarysize(const std::vector<double> &counts_hist);

double 
cl92_estimate_librarysize(const std::vector<double> &counts_hist);

/*
double
harris_lowerbound_librarysize(const std::vector<double> &counts_hist,
			      const double tolerance,
			      const size_t max_iter,
			      const size_t depth);
*/

double
my_harris(const bool VERBOSE,
	  const std::vector<double> &counts_hist,
	  const double tolerance,
	  const size_t max_iter,
	  const size_t depth);

/*double
upperbound_librarysize(const bool VERBOSE,
		       const std::vector<double> &counts_hist, 
		       size_t max_terms);
*/
#endif
