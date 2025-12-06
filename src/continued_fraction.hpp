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

#ifndef SRC_CONTINUED_FRACTION_HPP_
#define SRC_CONTINUED_FRACTION_HPP_

#include <cstddef>
#include <format>
#include <fstream>
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
  [[nodiscard]] double
  operator()(const double val) const;

  // Evaluate the continued fraction
  [[nodiscard]] double
  evaluate(const double val) const {
    return (*this)(val);
  };

  //////////////////////////////////////////
  // Extrapolation functions

  // Evaluate the continued fraction estimating distinct
  // along a curve from 0 to max_value
  void
  extrapolate_distinct(const double max_value, const double step_size,
                       std::vector<double> &estimates) const;

  bool
  is_valid() const {
    return !cf_coeffs.empty();
  }

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
  auto format(const ContinuedFraction &cf, auto &ctx) const {
    return std::formatter<std::string>::format(cf.tostring(), ctx);
  }
};

struct ContinuedFractionApproximation {
  // find best cont frac approx for estimating distinct
  [[nodiscard]] auto optimal_cf_distinct(
    const std::vector<double> &counts_hist) const -> ContinuedFraction;

  [[nodiscard]] auto get_diagonal() const -> int { return diagonal_idx; }

  int diagonal_idx{};       // the diagonal to work with for estimates
  std::size_t max_terms{};  // the maximum number of terms to try for a CF

  /* note: these never change */
  static constexpr std::size_t min_allowed_degree{4};

  // largest value to search for lowerbound and stability
  static constexpr double search_max_val{100.0};

  // step size for search of lowerbound and stability
  static constexpr double search_step_size{0.05};
};

[[nodiscard]] auto
check_yield_estimates_stability(const std::vector<double> &estimates) -> bool;

#endif  // SRC_CONTINUED_FRACTION_HPP_
