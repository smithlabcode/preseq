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

#include "continued_fraction.hpp"
#include "pade_approximant.hpp"

#include <vector>
#include <numeric>

double
chao87_lowerbound_librarysize(const std::vector<double> &counts_histogram);

double 
cl92_lowerbound_librarysize(const std::vector<double> &counts_histogram);

double
upperbound_librarysize(const std::vector<double> &counts_hist, 
		       size_t max_terms);

double
lowerbound_librarysize(const std::vector<double> &counts_hist,
		       const double upper_bound, //from upperbound_librarysize
                       const double step_size, const double max_val,
                       size_t max_terms)



#endif
