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

#include "continued_fraction.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

using std::fabs;
using std::isfinite;
using std::min;
using std::pow;
using std::vector;

// ADS: the std::pow function is used frequently to get (-1)^x for
// integer x. This doesn't make sense, and should be replaced at some
// point.

/* QUOTIENT DIFFERENCE ALGORITHM: compute continued fraction
 * coefficients vector for power series coefficients & vector for
 * continued fraction coefficients
 *
 * The negative sign for coefficients in the final loop is because we
 * evaluate a0/(1 + a1x/(1 + a2x/... while the algorithm is designed
 * for the a0/(1 - a1x/(1 - a2x/... see https://dlmf.nist.gov/3.10
 */
static void
quotdiff_algorithm(const vector<double> &ps_coeffs, vector<double> &cf_coeffs) {
  const size_t depth = ps_coeffs.size();  // degree of power series

  // q_table[0] never used, and undefined
  vector<vector<double>> q_table(depth, vector<double>(depth + 1, 0.0));
  // q_table[1][j] = ratio of ps coefficients
  for (size_t j = 0; j < depth - 1; j++)
    q_table[1][j] = ps_coeffs[j + 1] / ps_coeffs[j];

  // e_table[0] is always 0
  vector<vector<double>> e_table(depth, vector<double>(depth + 1, 0.0));
  // e_table[1] follows the general recurrence (same as in loop below)
  for (size_t j = 0; j < depth - 1; j++)
    e_table[1][j] = q_table[1][j + 1] - q_table[1][j] + e_table[0][j + 1];

  // using intial values of E(i)(j)'s and Q(i)(j)'s, fill rest of the
  // q table and e table
  for (size_t i = 2; i < depth; i++) {
    for (size_t j = 0; j < depth; j++)
      q_table[i][j] =
        q_table[i - 1][j + 1] * e_table[i - 1][j + 1] / e_table[i - 1][j];

    for (size_t j = 0; j < depth; j++)
      e_table[i][j] = q_table[i][j + 1] - q_table[i][j] + e_table[i - 1][j + 1];
  }

  cf_coeffs.resize(depth);
  // first CT coefficient is first PS coefficient
  cf_coeffs[0] = ps_coeffs[0];
  // set remaining CF coefficients from e and q table values
  for (size_t i = 1; i < depth; ++i)
    cf_coeffs[i] = (i % 2 == 0) ? -e_table[i / 2][0] : -q_table[(i + 1) / 2][0];
}

/* compute CF coeffs when upper_offset > 0 above the diagonal; this
 * means degree of polynomial in numerator of Pade approximant is
 * greater than degree of polynomial in the denominator
 */
static void
quotdiff_above_diagonal(const vector<double> &ps_coeffs, const size_t offset,
                        vector<double> &cf_coeffs,
                        vector<double> &offset_coeffs) {
  // get the high order PS coeffs for approximation by CF
  vector<double> high_ps_coeffs(begin(ps_coeffs) + offset, end(ps_coeffs));

  // use QD algorithm to determine CF coefficients
  quotdiff_algorithm(high_ps_coeffs, cf_coeffs);

  // first "offset" coeffs are equal to PS coeffs
  offset_coeffs = ps_coeffs;
  offset_coeffs.resize(offset);
}

// calculate CF coeffs when lower_offset > 0
static void
quotdiff_below_diagonal(const vector<double> &ps_coeffs, const size_t offset,
                        vector<double> &cf_coeffs,
                        vector<double> &offset_coeffs) {
  // need to work with reciprocal series g = 1/f, then invert
  vector<double> recip_ps_coeffs(ps_coeffs.size());
  recip_ps_coeffs[0] = 1.0 / ps_coeffs[0];
  for (size_t i = 1; i < ps_coeffs.size(); ++i) {
    double x = 0.0;
    for (size_t j = 0; j < i; ++j)
      x += ps_coeffs[i - j] * recip_ps_coeffs[j];

    recip_ps_coeffs[i] = -x / ps_coeffs[0];
  }

  // qd to compute cf_coeffs using remaining coeffs
  vector<double> high_recip_ps_coeffs(begin(recip_ps_coeffs) + offset,
                                      end(recip_ps_coeffs));
  quotdiff_algorithm(high_recip_ps_coeffs, cf_coeffs);

  // set offset coeffs to 1st "offset" PS coeffs of 1/f (reciprocal)
  offset_coeffs = recip_ps_coeffs;
  offset_coeffs.resize(offset);
}

void
truncate_degree(const size_t n_terms, ContinuedFraction &the_cf) {
  if (the_cf.degree < n_terms) {
    the_cf = ContinuedFraction();
  }
  else {
    the_cf.ps_coeffs.resize(n_terms);
    the_cf.cf_coeffs.resize(n_terms - the_cf.offset_coeffs.size());
    the_cf.degree = n_terms;
  }
}

ContinuedFraction::ContinuedFraction(const vector<double> &ps_cf, const int di,
                                     const size_t dg) :
  ps_coeffs(ps_cf), diagonal_idx(di), degree(dg) {
  if (diagonal_idx == 0)
    quotdiff_algorithm(ps_coeffs, cf_coeffs);
  else if (diagonal_idx > 0)
    quotdiff_above_diagonal(ps_coeffs, diagonal_idx, cf_coeffs, offset_coeffs);
  else  // if (cont_frac_estimate.lower_offset > 0) {
    quotdiff_below_diagonal(ps_coeffs, -diagonal_idx, cf_coeffs, offset_coeffs);
  // NOTE: negative sign "-" (-diagonal_idx > 0) for below diagonal
}

////////////////////////////////////////////////////////////////////////
//// FUNCTIONS TO EVALUATE CONTINUED FRACTIONS AT A POINT

static double
get_rescale_value(const double numerator, const double denominator) {
  static const double tolerance = 1e-20;  // magic
  const double rescale_val = fabs(numerator) + fabs(denominator);
  if (rescale_val > 1.0 / tolerance)
    return 1.0 / rescale_val;
  else if (rescale_val < tolerance)
    return 1.0 / rescale_val;
  return 1.0;
}

/* calculate ContinuedFraction approx when there is no offset uses euler's
 * recursion
 */
static double
evaluate_on_diagonal(const vector<double> &cf_coeffs, const double val,
                     const size_t depth) {
  // initialize
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;

  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0;

  for (size_t i = 1; i < min(cf_coeffs.size(), depth); i++) {
    // calculate current values
    current_num = prev_num1 + cf_coeffs[i] * val * prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i] * val * prev_denom2;

    // update previous values
    prev_num2 = prev_num1;
    prev_num1 = current_num;

    prev_denom2 = prev_denom1;
    prev_denom1 = current_denom;

    // now rescale all values
    const double rescale_val = get_rescale_value(current_num, current_denom);

    current_num *= rescale_val;
    current_denom *= rescale_val;

    prev_num1 *= rescale_val;
    prev_num2 *= rescale_val;

    prev_denom1 *= rescale_val;
    prev_denom2 *= rescale_val;
  }
  return current_num / current_denom;
}

static double
evaluate_power_series(const vector<double> &ps_coeffs, const double val) {
  double x = 0.0;
  for (size_t i = 0; i < ps_coeffs.size(); i++)
    x += ps_coeffs[i] * pow(val, i);
  return x;
}

/* evaluate CF when upper_offset > 0 using Euler's recursion */
static double
evaluate_above_diagonal(const vector<double> &cf_coeffs,
                        const vector<double> &offset_coeffs, const double val,
                        const size_t depth) {
  const double cf_part =
    evaluate_on_diagonal(cf_coeffs, val, depth - offset_coeffs.size());

  const double ps_part = evaluate_power_series(offset_coeffs, val);

  return ps_part + pow(val, offset_coeffs.size()) * cf_part;
}

// calculate ContinuedFraction approx when lower_offdiag > 0
static double
evaluate_below_diagonal(const vector<double> &cf_coeffs,
                        const vector<double> &offset_coeffs, const double val,
                        const size_t depth) {
  const double cf_part =
    evaluate_on_diagonal(cf_coeffs, val, depth - offset_coeffs.size());

  const double ps_part = evaluate_power_series(offset_coeffs, val);

  // recall that if lower_offset > 0, we are working with 1/f, invert approx
  return 1.0 / (ps_part + pow(val, offset_coeffs.size()) * cf_part);
}

// evaluate CF at a given point
double
ContinuedFraction::operator()(const double val) const {
  if (diagonal_idx > 0)
    return evaluate_above_diagonal(cf_coeffs, offset_coeffs, val, degree);
  else if (diagonal_idx < 0)
    return evaluate_below_diagonal(cf_coeffs, offset_coeffs, val, degree);
  else
    return evaluate_on_diagonal(cf_coeffs, val, degree);
}

std::ostream &
operator<<(std::ostream &the_stream, const ContinuedFraction &cf) {
  using std::ios_base;
  using std::setw;

  ios_base::fmtflags orig_flags = the_stream.flags();
  the_stream.setf(ios_base::fixed, ios_base::floatfield);
  the_stream.precision(2);
  the_stream << "OFFSET_COEFFS" << '\t' << "PS_COEFFS" << '\n';
  const size_t offset = cf.offset_coeffs.size();
  for (size_t i = 0; i < offset; ++i)
    the_stream << setw(12) << cf.offset_coeffs[i] << '\t' << setw(12)
               << cf.ps_coeffs[i] << '\n';
  the_stream << "CF_COEFFS" << '\n';
  for (size_t i = 0; i < cf.cf_coeffs.size(); ++i)
    the_stream << setw(12) << cf.cf_coeffs[i] << '\t' << setw(12)
               << cf.ps_coeffs[i + offset] << '\n';
  the_stream.flags(orig_flags);
  return the_stream;
}

// estimate yields by evaluating the CF at given points
void
ContinuedFraction::extrapolate_distinct(const double max_value,
                                        const double step_size,
                                        vector<double> &estimates) const {
  estimates.clear();
  estimates.push_back(0);
  for (double t = step_size; t <= max_value; t += step_size)
    estimates.push_back(t * operator()(t));
}

////////////////////////////////////////////////////////////////////////
////////////////  CONTINUED FRACTION APPROXIMATION CLASS BELOW

typedef ContinuedFractionApproximation CFA;

const size_t CFA::min_allowed_degree = 4;
const double CFA::search_max_val = 100;
const double CFA::search_step_size = 0.05;

/* check if a sequence of estimates are "stable": in [0, infty,
 * increasing, negative 2nd deriv
 */
bool
check_yield_estimates_stability(const vector<double> &estimates) {
  // require estimates are non-negative and finite
  for (size_t i = 0; i < estimates.size(); ++i)
    if (!std::isfinite(estimates[i]) || estimates[i] < 0.0)
      return false;

  // require estimate to be increasing
  for (size_t i = 1; i < estimates.size(); ++i)
    if (estimates[i] < estimates[i - 1])
      return false;

  // require negative second derivative
  for (size_t i = 2; i < estimates.size(); ++i)
    if (estimates[i - 1] - estimates[i - 2] < estimates[i] - estimates[i - 1])
      return false;

  return !estimates.empty();
}

/* Finds the optimal number of terms (i.e. degree, depth, etc.) of the
 * continued fraction by checking for stability of estimates at
 * specific points for yield. New way for searching for optimal CF
 */
ContinuedFraction
CFA::optimal_cont_frac_distinct(const vector<double> &counts_hist) const {
  // we expect to use an underestimate, but this is dealt with outside
  // by ensuring we have an even number of max terms

  if (max_terms >= counts_hist.size())
    return ContinuedFraction();

  vector<double> ps_coeffs;
  for (size_t j = 1; j <= max_terms; j++)
    ps_coeffs.push_back(counts_hist[j] * pow(-1.0, j + 1));

  ContinuedFraction full_cf(ps_coeffs, diagonal_idx, max_terms);

  // if max terms in {3,4,5,6}, check only that degree
  if (max_terms >= 3 && max_terms <= 6) {
    vector<double> estimates;
    full_cf.extrapolate_distinct(search_max_val, search_step_size, estimates);
    if (check_yield_estimates_stability(estimates))
      return full_cf;
  }
  else {
    // if max terms >= 7, start at 7 and check increasing cont frac's
    for (size_t i = 7 + (max_terms % 2 == 0); i <= max_terms; i += 2) {
      ContinuedFraction trunc_cf(full_cf);
      truncate_degree(i, trunc_cf);
      vector<double> estimates;
      trunc_cf.extrapolate_distinct(search_max_val, search_step_size,
                                    estimates);
      if (check_yield_estimates_stability(estimates))
        return trunc_cf;
    }
  }
  // no stable continued fraction: return null
  return ContinuedFraction();
}
