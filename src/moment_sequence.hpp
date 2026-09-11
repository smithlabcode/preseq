/* Copyright (C) 2013-2026 Andrew D. Smith and Timothy Daley
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_MOMENT_SEQUENCE_HPP_
#define SRC_MOMENT_SEQUENCE_HPP_

#include "nlohmann/json.hpp"

#include <cstddef>
#include <vector>

// test Hankel moment matrix to ensure the moment sequence is positive definite
auto
ensure_positive_definite_moment_sequence(std::vector<double> &moments,
                                         const double tolerance) -> std::size_t;

class MomentSequence {
public:
  MomentSequence() = default;
  explicit MomentSequence(const std::vector<double> &observed_moments);
  MomentSequence(const std::vector<double> &alpha,
                 const std::vector<double> &beta) : alpha{alpha}, beta{beta} {};

  // quadrature rules using QR on Jacobi matrix
  [[nodiscard]] auto
  lower_quadrature_rules(const std::size_t n_points, const double tolerance,
                         const std::size_t max_iter)
    -> std::tuple<std::vector<double>, std::vector<double>>;

private:
  // Estimate 3-term recurrence
  // these will be removed from the header when they are tested
  void
  unmodified_chebyshev();

  void
  full_three_term_recurrence(std::vector<double> &full_alpha,
                             std::vector<double> &full_beta);

  std::vector<double> moments;
  // 3-term recurrence
  std::vector<double> alpha;
  std::vector<double> beta;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(MomentSequence, moments, alpha, beta)
};

#endif  // SRC_MOMENT_SEQUENCE_HPP_
