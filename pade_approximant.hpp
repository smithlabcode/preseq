/*    Copyright (C) 2011 University of Southern California and
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

#ifndef PADE_APPROXIMANT_HPP
#define PADE_APPROXIMANT_HPP

#include <vector>

void 
compute_pade_coeffs(const std::vector<double> &coeffs, 
		    const size_t numer_size, const size_t denom_size,
		    std::vector<double> &num_coeffs, std::vector<double> &denom_coeffs);

double
compute_pade_approx_numerator(const double time,
			      const std::vector<double> &numerator);

double
compute_pade_approx_denominator(const double time,
				const std::vector<double> &denominator);

double
locate_polynomial_zero(const std::vector<double> &coeffs,
		       const double lower_lim, const double upper_lim,
		       const double tolerance);

void
compute_pade_curve(const std::vector<double> &coeffs,
		   const double max_time,
		   const double time_step,
		   const double allowable_defect_error,
		   const double tolerance,
		   const size_t denom_size,
		   const bool VERBOSE,
		   std::vector<double> &numerator_approx,
		   std::vector<double> &denominator_approx,
		   bool &defect_flag);

#endif
