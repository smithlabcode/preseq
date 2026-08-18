/* Copyright (C) 2013-2026 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Authors: Andrew D. Smith and Timothy Daley
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

#ifndef SRC_CONTINUED_FRACTION_HPP_
#define SRC_CONTINUED_FRACTION_HPP_

#include <cstddef>
#include <format>
#include <fstream>
#include <string>
#include <vector>

struct ContinuedFraction {
  // Constructors
  ContinuedFraction() = default;
  ContinuedFraction(const std::vector<double> &ps_cf, const int di,
                    const std::size_t dg);
  // Assumes diagonal is 0
  ContinuedFraction(const std::vector<double> &hist,
                    const std::size_t max_terms);

  // Evaluate the continued fraction
  [[nodiscard]] auto
  evaluate(const double val) const -> double;

  // Evaluate the continued fraction
  [[nodiscard]] auto
  operator()(const double val) const -> double {
    return evaluate(val);
  }

  // Extrapolation functions

  // Evaluate the continued fraction estimating distinct along a curve from 0
  // to max_value
  void
  extrapolate_distinct(const double max_value, const double step_size,
                       std::vector<double> &estimates) const;

  [[nodiscard]] auto
  is_valid() const -> bool {
    return !cf_coeffs.empty();
  }

  [[nodiscard]] auto
  return_degree() const -> std::size_t {
    return degree;
  }

  void
  extrapolate_curve(const double initial_distinct, const double vals_sum,
                    const double initial_sample_size, const double step_size,
                    const double max_sample_size,
                    std::vector<double> &estimates) const;

  [[nodiscard]] auto
  tostring() const -> std::string;

  std::vector<double> ps_coeffs;
  std::vector<double> cf_coeffs;
  std::vector<double> offset_coeffs;
  int diagonal_idx{};
  std::size_t degree{};
};

// get continued fraction with lower degree
void
decrease_degree(const std::size_t decrement, ContinuedFraction &cf);

void
truncate_degree(const std::size_t truncated_degree, ContinuedFraction &cf);

inline auto
operator<<(std::ostream &out, const ContinuedFraction &cf) -> std::ostream & {
  return out << cf.tostring();
}

template <>
struct std::formatter<ContinuedFraction> : std::formatter<std::string> {
  auto
  format(const ContinuedFraction &cf, auto &ctx) const {
    return std::formatter<std::string>::format(cf.tostring(), ctx);
  }
};

class ContinuedFractionApproximation {
public:
  ContinuedFractionApproximation(const int di, const std::size_t mt) :
    diagonal_idx{di}, max_terms{mt} {}

  // find best cont frac approx for estimating distinct
  [[nodiscard]] auto
  optimal_cf_distinct(const std::vector<double> &counts_hist) const
    -> ContinuedFraction;

  [[nodiscard]] auto
  get_diagonal() const -> int {
    return diagonal_idx;
  }

private:
  int diagonal_idx{};       // the diagonal to work with for estimates
  std::size_t max_terms{};  // the maximum number of terms to try for a CF

  /* note: these never change */
  static const std::size_t min_allowed_degree;

  // largest value to search for lowerbound and stability
  static const double search_max_val;

  // step size for search of lowerbound and stability
  static const double search_step_size;
};

[[nodiscard]] auto
check_yield_estimates_stability(const std::vector<double> &estimates) -> bool;

#endif  // SRC_CONTINUED_FRACTION_HPP_
