/* Copyright (C) 2013 University of Southern California and
 *                    Andrew D. Smith and Timothy Daley
 *
 * Authors: Andrew D. Smith and Timothy Daley
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_MOMENT_SEQUENCE_HPP_
#define SRC_MOMENT_SEQUENCE_HPP_

#include <cstddef>
#include <vector>

// test Hankel moment matrix to ensure the moment sequence
// is positive definite
std::size_t
ensure_pos_def_mom_seq(std::vector<double> &moments, const double tolerance);

struct MomentSequence {
  MomentSequence() {}
  explicit MomentSequence(const std::vector<double> &obs_moms);

  MomentSequence(const std::vector<double> &alpha,
                 const std::vector<double> &beta) : alpha{alpha}, beta{beta} {}

  // Estimate 3-term recurrence
  // these will be removed from the header when they are tested
  void
  unmodified_Chebyshev();

  void
  full_3term_recurrence(std::vector<double> &full_alpha,
                        std::vector<double> &full_beta);

  // quadrature rules using QR on Jacobi matrix
  bool
  lower_quadrature_rules(const std::size_t n_points, const double tolerance,
                         const std::size_t max_iter,
                         std::vector<double> &points,
                         std::vector<double> &weights);

  std::vector<double> moments;
  // 3-term recurrence
  std::vector<double> alpha;
  std::vector<double> beta;
};

#endif  // SRC_MOMENT_SEQUENCE_HPP_
