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

#ifndef SRC_CONTINUED_FRACTION_HPP_
#define SRC_CONTINUED_FRACTION_HPP_

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>

struct ContinuedFraction {
  // Constructors
  ContinuedFraction() : diagonal_idx(0), degree(0ul) {}
  ContinuedFraction(const std::vector<double> &ps_cf, const int di,
                    const size_t dg);

  // Evaluate the continued fraction
  double operator()(const double val) const;

  //////////////////////////////////////////
  // Extrapolation functions

  // Evaluate the continued fraction estimating distinct
  // along a curve from 0 to max_value
  void extrapolate_distinct(const double max_value, const double step_size,
                            std::vector<double> &estimates) const;

  bool is_valid() const { return !cf_coeffs.empty(); }
  size_t return_degree() const { return degree; }

  std::vector<double> ps_coeffs;
  std::vector<double> cf_coeffs;
  std::vector<double> offset_coeffs;
  int diagonal_idx;
  size_t degree;
};

// get continued fraction with lower degree
void
decrease_degree(const size_t decrement, ContinuedFraction &cf);
void
truncate_degree(const size_t truncated_degree, ContinuedFraction &cf);

std::ostream &
operator<<(std::ostream &out, const ContinuedFraction &cf);

class ContinuedFractionApproximation {
public:
  ContinuedFractionApproximation(const int di, const size_t mt) :
    diagonal_idx(di), max_terms(mt) {}

  // find best cont frac approx for estimating distinct
  ContinuedFraction
  optimal_cont_frac_distinct(const std::vector<double> &counts_hist) const;

  int get_diagonal() const { return diagonal_idx; }

private:
  int diagonal_idx;  // the diagonal to work with for estimates
  size_t max_terms;  // the maximum number of terms to try for a CF

  /* note: these never change */
  static const size_t min_allowed_degree;

  // largest value to search for lowerbound and stability
  static const double search_max_val;

  // step size for search of lowerbound and stability
  static const double search_step_size;
};

bool
check_yield_estimates_stability(const std::vector<double> &estimates);

#endif  // SRC_CONTINUED_FRACTION_HPP_
