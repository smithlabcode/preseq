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

#include "continued_fraction.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-avoid-magic-numbers,*-narrowing-conversions)

// ADS: the std::pow function is used frequently to get (-1)^x for integer
// x. This doesn't make sense, and should be replaced at some point.

/* QUOTIENT DIFFERENCE ALGORITHM: compute continued fraction
 * coefficients vector for power series coefficients & vector for
 * continued fraction coefficients
 *
 * The negative sign for coefficients in the final loop is because we
 * evaluate a0/(1 + a1x/(1 + a2x/... while the algorithm is designed
 * for the a0/(1 - a1x/(1 - a2x/... see https://dlmf.nist.gov/3.10
 */
[[nodiscard]] static auto
quotdiff_algorithm(const std::vector<double> &ps_coeffs)
  -> std::vector<double> {
  const std::size_t depth = std::size(ps_coeffs);  // degree of power series
  assert(depth > 0LU);

  // q_table[0] never used, and undefined
  auto q_table = std::vector(depth, std::vector<double>(depth + 1, 0.0));
  // q_table[1][j]: ratio of ps coefficients
  for (std::size_t j = 0; j + 1 < depth; ++j)
    q_table[1][j] = ps_coeffs[j + 1] / ps_coeffs[j];

  // e_table[0] is always 0
  auto e_table = std::vector(depth, std::vector<double>(depth + 1, 0.0));
  // e_table[1] follows the general recurrence (same as in loop below)
  for (std::size_t j = 0; j + 1 < depth; ++j)
    e_table[1][j] = q_table[1][j + 1] - q_table[1][j] + e_table[0][j + 1];

  // using intial values of E(i)(j)'s and Q(i)(j)'s, fill rest of the
  // q table and e table
  for (std::size_t i = 2; i < depth; ++i) {
    for (std::size_t j = 0; j < depth; ++j)
      q_table[i][j] =
        q_table[i - 1][j + 1] * e_table[i - 1][j + 1] / e_table[i - 1][j];

    for (std::size_t j = 0; j < depth; ++j)
      e_table[i][j] = q_table[i][j + 1] - q_table[i][j] + e_table[i - 1][j + 1];
  }

  std::vector<double> cf_coeffs(depth);
  // first CT coefficient is first PS coefficient
  cf_coeffs[0] = ps_coeffs[0];
  // set remaining CF coefficients from e and q table values
  for (std::size_t i = 1; i < depth; ++i)
    cf_coeffs[i] = (i % 2 == 0) ? -e_table[i / 2][0] : -q_table[(i + 1) / 2][0];
  return cf_coeffs;
}

/* compute CF coeffs when upper_offset > 0 above the diagonal; this
 * means degree of polynomial in numerator of Pade approximant is
 * greater than degree of polynomial in the denominator
 */
[[nodiscard]] static auto
quotdiff_above_diagonal(const std::vector<double> &ps_coeffs,
                        const std::size_t offset)
  -> std::tuple<std::vector<double>, std::vector<double>> {
  // get the high order PS coeffs for approximation by CF
  std::vector<double> high_ps_coeffs(std::cbegin(ps_coeffs) + offset,
                                     std::cend(ps_coeffs));
  // use QD algorithm to determine CF coefficients
  auto cf_coeffs = quotdiff_algorithm(high_ps_coeffs);
  // first "offset" coeffs are equal to PS coeffs
  auto offset_coeffs = ps_coeffs;
  offset_coeffs.resize(offset);
  return std::tuple{std::move(cf_coeffs), std::move(offset_coeffs)};
}

// calculate CF coeffs when lower_offset > 0
[[nodiscard]] static auto
quotdiff_below_diagonal(const std::vector<double> &ps_coeffs,
                        const std::size_t offset)
  -> std::tuple<std::vector<double>, std::vector<double>> {
  // need to work with reciprocal series g = 1/f, then invert
  std::vector<double> recip_ps_coeffs(std::size(ps_coeffs));
  recip_ps_coeffs[0] = 1.0 / ps_coeffs[0];
  for (std::size_t i = 1; i < std::size(ps_coeffs); ++i) {
    double x = 0.0;
    for (std::size_t j = 0; j < i; ++j)
      x += ps_coeffs[i - j] * recip_ps_coeffs[j];

    recip_ps_coeffs[i] = -x / ps_coeffs[0];
  }

  // qd to compute cf_coeffs using remaining coeffs
  std::vector<double> high_recip_ps_coeffs(
    std::cbegin(recip_ps_coeffs) + offset, std::cend(recip_ps_coeffs));
  auto cf_coeffs = quotdiff_algorithm(high_recip_ps_coeffs);

  // set offset coeffs to 1st "offset" PS coeffs of 1/f (reciprocal)
  auto offset_coeffs = recip_ps_coeffs;
  offset_coeffs.resize(offset);
  return std::tuple{std::move(cf_coeffs), std::move(offset_coeffs)};
}

/* decrease degree of CF keeping coeffs equal to original */
void
decrease_degree(const std::size_t decrement, ContinuedFraction &cf) {
  assert(decrement < cf.degree);
  cf.ps_coeffs.resize(std::size(cf.ps_coeffs) - decrement);
  cf.cf_coeffs.resize(std::size(cf.cf_coeffs) - decrement);
  cf.degree -= decrement;
}

void
truncate_degree(const std::size_t n_terms, ContinuedFraction &cf) {
  if (cf.degree < n_terms) {
    cf = ContinuedFraction();
  }
  else {
    cf.ps_coeffs.resize(n_terms);
    cf.cf_coeffs.resize(n_terms - std::size(cf.offset_coeffs));
    cf.degree = n_terms;
  }
}

ContinuedFraction::ContinuedFraction(const std::vector<double> &ps_cf,
                                     const int di, const std::size_t dg) :
  ps_coeffs(ps_cf),
  diagonal_idx(di), degree(dg) {
  if (diagonal_idx == 0)
    cf_coeffs = quotdiff_algorithm(ps_coeffs);
  else if (diagonal_idx > 0)
    std::tie(cf_coeffs, offset_coeffs) =
      quotdiff_above_diagonal(ps_coeffs, diagonal_idx);
  else  // if (cont_frac_estimate.lower_offset > 0) {
    std::tie(cf_coeffs, offset_coeffs) =
      quotdiff_below_diagonal(ps_coeffs, -diagonal_idx);
  // NOTE: negative sign "-" (-diagonal_idx > 0) for below diagonal
}

ContinuedFraction::ContinuedFraction(const std::vector<double> &hist,
                                     const std::size_t max_terms) :
  degree{max_terms} {
  for (std::size_t j = 1; j <= max_terms; ++j)
    ps_coeffs.push_back(hist[j] * std::pow(-1.0, j + 1));
  cf_coeffs = quotdiff_algorithm(ps_coeffs);
}

/// Functions to evaluate continued fractions at a point

[[nodiscard]] static auto
get_rescale_value(const double numerator, const double denominator) -> double {
  static const double tolerance = 1e-20;  // magic
  const double rescale_val = std::fabs(numerator) + std::fabs(denominator);
  if (rescale_val > 1.0 / tolerance)
    return 1.0 / rescale_val;
  else if (rescale_val < tolerance)
    return 1.0 / rescale_val;
  return 1.0;
}

/* calculate ContinuedFraction approx when there is no offset uses euler's
 * recursion
 */
[[nodiscard]] static auto
evaluate_on_diagonal(const std::vector<double> &cf_coeffs, const double val,
                     const std::size_t depth) -> double {
  // initialize
  double current_numer{};
  double prev_numer1 = cf_coeffs[0];
  double prev_numer2{};

  double current_denom{};
  double prev_denom1{1.0};
  double prev_denom2{1.0};

  const auto lim = std::min(std::size(cf_coeffs), depth);
  for (std::size_t i = 1; i < lim; ++i) {
    // calculate current values
    current_numer = prev_numer1 + cf_coeffs[i] * val * prev_numer2;
    current_denom = prev_denom1 + cf_coeffs[i] * val * prev_denom2;

    // update previous values
    prev_numer2 = prev_numer1;
    prev_numer1 = current_numer;

    prev_denom2 = prev_denom1;
    prev_denom1 = current_denom;

    // now rescale all values
    const double rescale_val = get_rescale_value(current_numer, current_denom);

    current_numer *= rescale_val;
    current_denom *= rescale_val;

    prev_numer1 *= rescale_val;
    prev_numer2 *= rescale_val;

    prev_denom1 *= rescale_val;
    prev_denom2 *= rescale_val;
  }
  return current_numer / current_denom;
}

[[nodiscard]] static auto
evaluate_power_series(const std::vector<double> &ps_coeffs, const double val)
  -> double {
  double x{};
  for (std::size_t i = 0; i < std::size(ps_coeffs); ++i)
    x += ps_coeffs[i] * std::pow(val, i);
  return x;
}

/* evaluate CF when upper_offset > 0 using Euler's recursion */
[[nodiscard]] static auto
evaluate_above_diagonal(const std::vector<double> &cf_coeffs,
                        const std::vector<double> &offset_coeffs,
                        const double val, const std::size_t depth) -> double {
  const double cf_part =
    evaluate_on_diagonal(cf_coeffs, val, depth - std::size(offset_coeffs));

  const double ps_part = evaluate_power_series(offset_coeffs, val);

  return ps_part + std::pow(val, std::size(offset_coeffs)) * cf_part;
}

// calculate ContinuedFraction approx when lower_offdiag > 0
[[nodiscard]] static auto
evaluate_below_diagonal(const std::vector<double> &cf_coeffs,
                        const std::vector<double> &offset_coeffs,
                        const double val, const std::size_t depth) -> double {
  const double cf_part =
    evaluate_on_diagonal(cf_coeffs, val, depth - std::size(offset_coeffs));

  const double ps_part = evaluate_power_series(offset_coeffs, val);

  // recall that if lower_offset > 0, we are working with 1/f, invert approx
  return 1.0 / (ps_part + std::pow(val, std::size(offset_coeffs)) * cf_part);
}

// evaluate CF at a given point
[[nodiscard]] auto
ContinuedFraction::evaluate(const double val) const -> double {
  if (diagonal_idx > 0)
    return evaluate_above_diagonal(cf_coeffs, offset_coeffs, val, degree);
  else if (diagonal_idx < 0)
    return evaluate_below_diagonal(cf_coeffs, offset_coeffs, val, degree);
  else
    return evaluate_on_diagonal(cf_coeffs, val, degree);
}

void
ContinuedFraction::extrapolate_curve(const double initial_distinct,
                                     const double n_observations,
                                     const double initial_sample_size,
                                     const double step_size,
                                     const double max_sample_size,
                                     std::vector<double> &estimates) const {
  double current_sample_size = initial_sample_size;
  while (current_sample_size < max_sample_size) {
    const double fold = (current_sample_size - n_observations) / n_observations;
    assert(fold >= 0.0);
    estimates.push_back(initial_distinct + fold * evaluate(fold));
    current_sample_size += step_size;
  }
}

[[nodiscard]] auto
ContinuedFraction::tostring() const -> std::string {
  std::ostringstream the_stream;
  std::ios_base::fmtflags orig_flags = the_stream.flags();
  the_stream.setf(std::ios_base::fixed, std::ios_base::floatfield);
  the_stream.precision(2);
  the_stream << "OFFSET_COEFFS" << '\t' << "PS_COEFFS" << '\n';
  const std::size_t offset = std::size(offset_coeffs);
  for (std::size_t i = 0; i < offset; ++i)
    the_stream << std::setw(12) << offset_coeffs[i] << '\t' << std::setw(12)
               << ps_coeffs[i] << '\n';
  the_stream << "CF_COEFFS" << '\n';
  for (std::size_t i = 0; i < std::size(cf_coeffs); ++i)
    the_stream << std::setw(12) << cf_coeffs[i] << '\t' << std::setw(12)
               << ps_coeffs[i + offset] << '\n';
  the_stream.flags(orig_flags);
  return the_stream.str();
}

// estimate yields by evaluating the CF at given points
[[nodiscard]] auto
ContinuedFraction::extrapolate_distinct(const double max_value,
                                        const double step_size) const
  -> std::vector<double> {
  std::vector<double> estimates;
  estimates.push_back(0);
  auto t = step_size;
  while (t <= max_value) {
    estimates.push_back(t * evaluate(t));
    t += step_size;
  }
  return estimates;
}

/// Continued fraction _approximation_

/* check if a sequence of estimates are "stable": in [0, infty), increasing,
 * negative 2nd deriv
 */
[[nodiscard]] auto
check_yield_estimates_stability(const std::vector<double> &estimates) -> bool {
  // require estimates are non-negative and finite
  for (const auto estimate : estimates)
    if (!std::isfinite(estimate) || estimate < 0.0)
      return false;

  // require estimate to be increasing
  for (std::size_t i = 1; i < std::size(estimates); ++i)
    if (estimates[i] < estimates[i - 1])
      return false;

  // require negative second derivative
  const auto neg_2nd_deriv = [](const auto x1, const auto x2, const auto x3) {
    return (x2 - x1) < (x3 - x2);
  };
  for (std::size_t i = 2; i < std::size(estimates); ++i)
    if (neg_2nd_deriv(estimates[i - 2], estimates[i - 1], estimates[i]))
      return false;

  return !estimates.empty();
}

/* Finds the optimal number of terms (i.e. degree, depth, etc.) of the continued
 * fraction by checking for stability of estimates at specific points for
 * yield. New way for searching for optimal CF.
 */
[[nodiscard]] auto
ContinuedFractionApproximation::optimal_cf_distinct(
  const std::vector<double> &counts_hist) const -> ContinuedFraction {
  // we expect to use an underestimate, but this is dealt with outside
  // by ensuring we have an even number of max terms

  if (max_terms >= std::size(counts_hist))
    return {};

  std::vector<double> ps_coeffs;
  for (std::size_t j = 1; j <= max_terms; ++j)
    ps_coeffs.push_back(counts_hist[j] * std::pow(-1.0, j + 1));

  ContinuedFraction full_cf(ps_coeffs, diagonal_idx, max_terms);

  // if max terms in {3,4,5,6}, check only that degree
  if (max_terms >= 3 && max_terms <= 6) {
    const auto estimates =
      full_cf.extrapolate_distinct(search_max_val, search_step_size);
    if (check_yield_estimates_stability(estimates))
      return full_cf;
  }
  else {
    // if max terms >= 7, start at 7 and check increasing cont frac's
    for (std::size_t i = 7 + (max_terms % 2 == 0); i <= max_terms; i += 2) {
      ContinuedFraction truncated_cf(full_cf);
      truncate_degree(i, truncated_cf);
      const auto estimates =
        truncated_cf.extrapolate_distinct(search_max_val, search_step_size);
      if (check_yield_estimates_stability(estimates))
        return truncated_cf;
    }
  }
  // no stable continued fraction: return null
  return {};
}

// NOLINTEND(*-avoid-magic-numbers,*-narrowing-conversions)
