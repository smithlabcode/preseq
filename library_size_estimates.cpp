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

#include "library_size_estimates.hpp"

#include "pade_approximant.hpp"
#include "continued_fraction.hpp"
#include "smithlab_utils.hpp"

using std::string;
using std::vector;
using std::max;

using smithlab::log_sum_log_vec;

static const size_t MIN_ALLOWED_DEGREE = 6;

// compute the upperbound on library size by noting that p(x)/q(x)
// is a constant in the limit if p & q are of equal degree.
// furthermore if p & q are equal degree then the library_yield
// estimates are liberal, i.e. they overestimate the library_yield
// on average
double
upperbound_librarysize(const vector<double> &counts_hist, size_t max_terms) {
  // need max_terms = L+M+1 to be even so that L+M is odd so that we
  //can take lim_{t \to \infty} [L+1, M]
  if (max_terms % 2 == 1)
    --max_terms;
  
  vector<double> coeffs;
  for(size_t i = 0; i < max_terms; i++)
    coeffs.push_back(pow(-1, i+2)*counts_hist[i+1]);

  bool ACCEPT_UPPER_BOUND = false;
  double upper_bound;

  while (max_terms >= MIN_ALLOWED_DEGREE
	 && !ACCEPT_UPPER_BOUND) {
    vector<double> denom_vec;
    vector<double> num_vec;
    const size_t numer_size = max_terms/2;  
    const size_t denom_size = max_terms - numer_size;
    // numer_size = degree(p)+1, denom_size = degree(q)
    // library_yield = x*p(x)/q(x), ensure that degree(x*p)=degree(q)
    assert(numer_size == denom_size);

    // consider upper bound if pade approx is acceptable
    const bool ACCEPT_APPROX = 
      compute_pade_coeffs(coeffs, numer_size, denom_size, num_vec, denom_vec); 

    // lim(xp(x)/q(x)) = p_{numer_size-1}/q_{denom_size}
    //coefficients are in order of degree
    upper_bound = num_vec.back()/denom_vec.back();
    ACCEPT_UPPER_BOUND = ACCEPT_APPROX && upper_bound > 0;
      max_terms -= 2;
  }
  return upper_bound;
}

// library_yield = xp(x)/q(x),
// so if degree(q) > degree(p)+1, then library_yield acts 
//like 1/x^n for some n > 0 in the limit and therefore goes to zero
//since it approximates the library yield in the neighborhood of zero
//there is global max, so if we choose a conservative approx, this is 
//a lower bound on library_size
double
lowerbound_librarysize(const vector<double> &counts_hist,
		       const double upper_bound, //from upper_bound_librarysize(.)
                       const double step_size, const double max_val,
                       size_t max_terms) {
  // the derivative must always be less than the number of distinct reads
  //in the initial sample
  const double distinct_reads = 
    accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

  // make sure we are using appropriate order estimate
  if (max_terms % 2 == 0)
    --max_terms;
  
  // Initialize the set of coefficients for the original power series
  vector<double> coeffs(max_terms, 0.0);
  for (size_t j = 0; j < max_terms; ++j)
    coeffs[j] = counts_hist[j + 1]*pow(-1, j + 2);
  
  //iterate over max_terms to find largest local max as lower bound
  //theortically larger max_terms will be better approximations
  //==> larger lower bounds
  vector<double> possible_maxima_loc, possible_maxima;
  size_t n_terms = max_terms;
  while (n_terms > ContFracApprox::MINIMUM_ALLOWED_DEGREE) {
    
    cont_frac cf_estimate(coeffs, 2, 0);
    ContFracApprox CFestimator(cf_estimate, n_terms);
    // evaluate for the current "max_terms"
    possible_maxima.push_back(CFestimator.locate_local_max(0.0, max_val, 
							   step_size,
							   upper_bound, 
							   distinct_reads));

    possible_maxima_loc.push_back(CFestimator.evaluate(possible_maxima.back()));

    //move down in terms
    n_terms -= 2;
  }
  
  // Now compare the different degrees to see which is has the best
  // maxima
  double global_max = 0.0, global_max_loc = 0.0;
  for (size_t i = 0; i < possible_maxima_loc.size(); ++i)
    if (possible_maxima[i] > global_max && possible_maxima[i] < upper_bound) {
      global_max = possible_maxima[i];
      global_max_loc = possible_maxima_loc[i];
    }
  
  return global_max;
}

// Chao (Biometrics 1987) lower bound
double
chao87_lowerbound_librarysize(const vector<double> &counts_hist) {
  assert(counts_hist.size() >= 2);
  return accumulate(counts_hist.begin(), counts_hist.end(), 0.0) +
    counts_hist[1]*counts_hist[1]/(2.0*counts_hist[2]);
}


//Chao & Lee (JASA 1992) lower bound
double 
cl92_lowerbound_librarysize(const vector<double> &counts_hist) {
  
  double sample_size = 0.0;
  for(size_t i = 0; i <  counts_hist.size(); i++)
    sample_size += counts_hist[i]*i;
  const double distinct_vals = accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

  const double estim_coverage = 1.0 - counts_hist[1]/sample_size;
  const double naive_lowerbound = distinct_vals/estim_coverage;
  
  vector<double> log_coeff_var;
  for(size_t i = 2; i < counts_hist.size(); ++i)
    if (counts_hist[i] > 0)
      log_coeff_var.push_back(log(i) + log(i - 1) + log(counts_hist[i]));
  
  const double coeff_var = 
    naive_lowerbound*exp(log_sum_log_vec(log_coeff_var, log_coeff_var.size())
			 -log(sample_size) - log(sample_size - 1)) - 1;
  return(naive_lowerbound + counts_hist[1]*coeff_var/estim_coverage);
}


/* Chao & Lee (JASA 1992) lower bound
 */
// static double 
// compute_cl92_library_size_lower_bound(const vector<double> &vals_hist) {
  
//   double sample_size = 0.0;
//   for (size_t i = 0; i < vals_hist.size(); ++i)
//     sample_size += i*vals_hist[i]; 
//   const double distinct = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
  
//   // The "coverage" below is an estimate??
//   const double coverage = 1.0 - vals_hist[1]/sample_size;
//   const double naive_lowerbound = distinct/coverage;
  
//   /// This is sum(i(i-1)f_i) from 2.12 and 2.13 of CL92 paper
//   vector<double> log_cv_terms;
//   for (size_t i = 2; i < vals_hist.size(); ++i) 
//     if (vals_hist[i] > 0)
//       log_cv_terms.push_back((log(i) + log(i - 1)) + log(vals_hist[i]));

//   const double cv_terms_sum = (log_cv_terms.empty()) ? 0.0 : 
//     exp(log_sum_log_vec(log_cv_terms, log_cv_terms.size()));
  
//   const double coeff_variation_est =  max(naive_lowerbound*cv_terms_sum/
// 					  (sample_size*(sample_size - 1.0)) - 1.0, 0.0);
  
//   // Corrected coefficient of variation below is from 2.13 in Chao &
//   // Lee (JASA 1992)
//   const double cv_bias_correction = 1.0 + 
//     ((1.0 - coverage)*cv_terms_sum)/((sample_size - 1.0)*coverage);
  
//   const double coeff_variation_corr = 
//     max(coeff_variation_est*cv_bias_correction,	0.0);
  
//   return naive_lowerbound + 
//     coeff_variation_corr*(sample_size*(1.0 - coverage)/coverage);
// }

// //Chao & Lee (JASA 1992) lower bound
// static double 
// compute_cl92_library_size_lower_bound(const vector<double> &vals_hist) {
//   double sample_size = 0.0;
//   for(size_t i = 0; i < vals_hist.size(); i++)
//     sample_size += i*vals_hist[i]; 
//   const double distinct = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
//   double coverage = 1 - vals_hist[1]/sample_size;
//   double naive_lowerbound = distinct/coverage;
//   vector<double> log_cv_terms;
//   for (size_t i = 2; i < vals_hist.size(); ++i)
//     if (vals_hist[i] > 0)
//       log_cv_terms.push_back(log(naive_lowerbound) + log(i) + log(i-1)
//                              -log(sample_size) - log(sample_size-1));
  
//   double coeff_variation = 0.0;
//   if (log_cv_terms.size() > 0)
//     coeff_variation = max(exp(log_sum_log_vec(log_cv_terms, 
// 					      log_cv_terms.size())) - 1, 0.0);
  
//   for(size_t i = 0; i < log_cv_terms.size(); i++)
//     log_cv_terms[i] -= log(naive_lowerbound);
//   const double corrected_coeff_variation = 
//     coeff_variation*(1+sample_size*(1 - coverage)*
// 		     exp(log_sum_log_vec(log_cv_terms, log_cv_terms.size()))/coverage);
  
//   return naive_lowerbound + 
//     sample_size*(1 - coverage)*corrected_coeff_variation/coverage;
// }
