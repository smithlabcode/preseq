/*    Copyright (C) 2013 University of Southern California and
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

#ifndef MOMENT_SEQUENCE_HPP
#define MOMENT_SEQUENCE_HPP

#include <vector>
#include <numeric>

// test Hankel moment matrix to ensure the moment sequence
// is positive definite
size_t ensure_pos_def_mom_seq(std::vector<double> &moments,
			      const double tolerance,
			      const bool VERBOSE);

struct MomentSequence {

  // Constructors
  MomentSequence() {}
  MomentSequence(const std::vector<double> &obs_moms);

  MomentSequence(const std::vector<double> &a,
		 const std::vector<double> &b):
    alpha(a), beta(b) {};



  // Estimate 3-term recurrence
  // these will be removed from the header when they are tested
  void unmodified_Chebyshev(const bool VERBOSE);

  void full_3term_recurrence(const bool VERBOSE,
			     std::vector<double> &full_alpha,
			     std::vector<double> &full_beta);

  // quadrature rules using QR on Jacobi matrix 
  bool Lower_quadrature_rules(const bool VERBOSE,
				 const size_t n_points,
				 const double tolerance, 
				 const size_t max_iter,
				 std::vector<double> &points,
				 std::vector<double> &weights);


  std::vector<double> moments;
  // 3-term recurrence
  std::vector<double> alpha;
  std::vector<double> beta;
};


#endif
