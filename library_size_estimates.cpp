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

// computes bounds on library size
// including lower bounds of Chao(Biometrics 1987) and 
// Chao & Lee (JASA 1992)


#include <numeric>
#include <vector>
#include <string>
#include <iostream>
#include <ostream>
#include <cassert>

#include "smithlab_utils.hpp"

#include "library_size_estimates.hpp"
#include "pade_approximant.hpp"
#include "continued_fraction.hpp"

using std::string;
using std::vector;
using std::max;

using smithlab::log_sum_log_vec;

static const size_t MIN_ALLOWED_DEGREE = 6;

/* compute the upperbound on library size by noting that p(x)/q(x) is
 * a constant in the limit if p & q are of equal degree. furthermore
 * if p & q are equal degree then the library_yield estimates are
 * liberal, i.e. they overestimate the library_yield on average
 */
double
upperbound_librarysize(const vector<double> &counts_hist, size_t max_terms) {
  // need max_terms = L + M + 1 to be even so that L + M is odd so
  // that we can take lim_{t \to \infty} [L+1, M]
  
  if (max_terms % 2 == 1)
    --max_terms;
  
  vector<double> ps_coeffs;
  for (size_t i = 0; i < max_terms; ++i)
    ps_coeffs.push_back(pow(-1, i + 2)*counts_hist[i + 1]);
  
  for (; max_terms >= MIN_ALLOWED_DEGREE; max_terms -= 2) {
    const size_t numer_degree = max_terms/2;  
    const size_t denom_degree = max_terms - numer_degree;
    assert(numer_degree == denom_degree);
    // numer_degree = degree(p)+1, denom_degree = degree(q); and
    // library_yield = x*p(x)/q(x), ensure that degree(x*p)=degree(q)
    
    // consider upper bound if pade approx is acceptable
    vector<double> numers, denoms;
    const bool accept_approx = 
      compute_pade_coeffs(ps_coeffs, numer_degree, denom_degree, numers, denoms); 
    
    // lim(xp(x)/q(x)) = p_{numer_size-1}/q_{denom_size} coefficients
    // are in order of degree
    const double upper_bound = numers.back()/denoms.back();
    if (accept_approx && upper_bound > 0.0 && std::isfinite(upper_bound))
      return upper_bound;
  }
  return -std::numeric_limits<double>::max();
}


// Chao (Biometrics 1987) lower bound
double
chao87_lowerbound_librarysize(const vector<double> &counts_hist) {
  assert(counts_hist.size() >= 2);
  return accumulate(counts_hist.begin(), counts_hist.end(), 0.0) +
    counts_hist[1]*counts_hist[1]/(2.0*counts_hist[2]);
}


// Chao & Lee (JASA 1992) lower bound
double 
cl92_lowerbound_librarysize(const vector<double> &counts_hist) {
  
  double sample_size = 0.0;
  for(size_t i = 0; i <  counts_hist.size(); i++)
    sample_size += i*counts_hist[i];
  
  const double distinct_vals = 
    accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

  const double estim_coverage = 1.0 - counts_hist[1]/sample_size;
  const double naive_lowerbound = distinct_vals/estim_coverage;
  
  vector<double> log_coeff_var;
  for (size_t i = 2; i < counts_hist.size(); ++i)
    if (counts_hist[i] > 0)
      log_coeff_var.push_back(log(i) + log(i - 1) + log(counts_hist[i]));
  
  const double coeff_var = 
    naive_lowerbound*exp(log_sum_log_vec(log_coeff_var, log_coeff_var.size())
			 -log(sample_size) - log(sample_size - 1)) - 1;
  
  return naive_lowerbound + counts_hist[1]*coeff_var/estim_coverage;
}
