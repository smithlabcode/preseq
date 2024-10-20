/* Copyright (C) 2013-2024 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Authors: Timothy Daley and Andrew Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "common.hpp"

#include "continued_fraction.hpp"

#include <unistd.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using std::array;
using std::begin;
using std::cbegin;
using std::cend;
using std::cerr;
using std::end;
using std::endl;
using std::min;
using std::mt19937;
using std::runtime_error;
using std::size_t;
using std::string;
using std::uint32_t;
using std::vector;

double
GoodToulmin2xExtrap(const vector<double> &counts_hist) {
  double two_fold_extrap = 0.0;
  for (size_t i = 0; i < counts_hist.size(); i++)
    two_fold_extrap += pow(-1.0, i + 1) * counts_hist[i];
  return two_fold_extrap;
}

// Lanczos approximation for gamma function for x >= 0.5 - essentially an
// approximation for (x-1)!
double
factorial(double x) {
  // constants
  static constexpr double LogRootTwoPi = 0.9189385332046727;
  static constexpr double Euler = 2.71828182845904523536028747135;
  array<double, 9> Lanczos{0.99999999999980993227684700473478,
                           676.520368121885098567009190444019,
                           -1259.13921672240287047156078755283,
                           771.3234287776530788486528258894,
                           -176.61502916214059906584551354,
                           12.507343278686904814458936853,
                           -0.13857109526572011689554707,
                           9.984369578019570859563e-6,
                           1.50563273514931155834e-7};

  // Approximation for factorial is actually x-1
  x -= 1.0;

  double Ag = Lanczos[0];
  for (auto k = 1u; k < size(Lanczos); k++)
    Ag += Lanczos[k] / (x + k);

  const double term1 = (x + 0.5) * log((x + 7.5) / Euler);
  const double term2 = LogRootTwoPi + log(Ag);

  return term1 + (term2 - 7.0);
}

// interpolate by explicit calculating the expectation
// for sampling without replacement;
// see K.L Heck 1975
// N total sample size; S the total number of distincts
// n sub sample size
double
interpolate_distinct(const vector<double> &hist, const size_t N, const size_t S,
                     const size_t n) {
  const double denom =
    factorial(N + 1) - factorial(n + 1) - factorial(N - n + 1);

  vector<double> numer(hist.size(), 0);
  for (size_t i = 1; i < hist.size(); i++) {
    // N - i -n + 1 should be greater than 0
    if (N < i + n) {
      numer[i] = 0;
    }
    else {
      const double x =
        (factorial(N - i + 1) - factorial(n + 1) - factorial(N - i - n + 1));
      numer[i] = exp(x - denom) * hist[i];
    }
  }
  return S - accumulate(cbegin(numer), cend(numer), 0);
}

static void
extrapolate_curve(const ContinuedFraction &the_cf,
                  const double initial_distinct, const double vals_sum,
                  const double initial_sample_size, const double step_size,
                  const double max_sample_size, vector<double> &estimates) {
  double curr_samp_sz = initial_sample_size;
  while (curr_samp_sz < max_sample_size) {
    const double fold = (curr_samp_sz - vals_sum) / vals_sum;
    assert(fold >= 0.0);
    estimates.push_back(initial_distinct + fold * the_cf(fold));
    curr_samp_sz += step_size;
  }
}

bool
extrap_single_estimate(const bool VERBOSE, const bool allow_defects,
                       const vector<double> &hist, size_t max_terms,
                       const int diagonal, const double step_size,
                       const double max_extrap,
                       vector<double> &yield_estimate) {
  yield_estimate.clear();

  const double vals_sum = get_counts_from_hist(hist);
  const double initial_distinct = accumulate(cbegin(hist), cend(hist), 0.0);

  // interpolate complexity curve by random sampling w/out replacement
  const size_t upper_limit = vals_sum;
  const size_t step = step_size;
  size_t sample = static_cast<size_t>(step_size);
  for (; sample < upper_limit; sample += step)
    yield_estimate.push_back(
      interpolate_distinct(hist, upper_limit, initial_distinct, sample));

  // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
  size_t first_zero = 1;
  while (first_zero < hist.size() && hist[first_zero] > 0)
    ++first_zero;

  // Ensure we are not using a zero term
  max_terms = min(max_terms, first_zero - 1);

  // refit curve for lower bound (degree of approx is 1 less than
  // max_terms)
  max_terms = max_terms - (max_terms % 2 == 1);

  if (allow_defects) {
    vector<double> ps_coeffs;
    for (size_t j = 1; j <= max_terms; j++)
      ps_coeffs.push_back(hist[j] * std::pow(-1.0, j + 1));

    const ContinuedFraction defect_cf(ps_coeffs, diagonal, max_terms);

    extrapolate_curve(defect_cf, initial_distinct, vals_sum, sample, step_size,
                      max_extrap, yield_estimate);

    if (VERBOSE)
      cerr << defect_cf << endl;
    // NO FAIL! defect mode doesn't care about failure
  }
  else {
    const ContinuedFractionApproximation lower_cfa(diagonal, max_terms);
    const ContinuedFraction lower_cf(
      lower_cfa.optimal_cont_frac_distinct(hist));

    // extrapolate curve
    if (lower_cf.is_valid()) {
      extrapolate_curve(lower_cf, initial_distinct, vals_sum, sample, step_size,
                        max_extrap, yield_estimate);
    }
    else {
      // FAIL! lower_cf unacceptable, need to bootstrap to obtain
      // estimates
      return false;
    }

    if (VERBOSE)
      cerr << lower_cf << endl;
  }
  // SUCCESS!!
  return true;
}

void
extrap_bootstrap(const bool VERBOSE, const bool allow_defects,
                 const uint32_t seed, const vector<double> &orig_hist,
                 const size_t n_bootstraps, const size_t orig_max_terms,
                 const int diagonal, const double bin_step_size,
                 const double max_extrap, const size_t max_iter,
                 vector<vector<double>> &bootstrap_estimates) {
  // clear returning vectors
  bootstrap_estimates.clear();

  // setup rng
  mt19937 rng(seed);

  const double initial_distinct =
    std::accumulate(cbegin(orig_hist), cend(orig_hist), 0.0);

  vector<size_t> orig_hist_distinct_counts;
  vector<double> distinct_orig_hist;
  for (size_t i = 0; i < orig_hist.size(); i++)
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }

  for (size_t iter = 0;
       (iter < max_iter && bootstrap_estimates.size() < n_bootstraps); ++iter) {
    if (VERBOSE && iter > 0 && iter % 72 == 0)
      cerr << endl;  // bootstrap success progress only 72 char wide

    vector<double> yield_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);

    const double sample_vals_sum = get_counts_from_hist(hist);

    // resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    // compute complexity curve by random sampling w/out replacement
    const size_t distinct = accumulate(cbegin(hist), cend(hist), 0.0);
    size_t curr_sample_sz = bin_step_size;
    while (curr_sample_sz < sample_vals_sum) {
      yield_vector.push_back(
        interpolate_distinct(hist, sample_vals_sum, distinct, curr_sample_sz));
      curr_sample_sz += bin_step_size;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t first_zero = 1;
    while (first_zero < hist.size() && hist[first_zero] > 0)
      ++first_zero;

    size_t max_terms = min(orig_max_terms, first_zero - 1);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 1);

    bool successful_bootstrap = false;
    // defect mode, simple extrapolation
    if (allow_defects) {
      vector<double> ps_coeffs;
      for (size_t j = 1; j <= max_terms; j++)
        ps_coeffs.push_back(hist[j] * std::pow(-1.0, j + 1));

      const ContinuedFraction defect_cf(ps_coeffs, diagonal, max_terms);

      extrapolate_curve(defect_cf, initial_distinct, sample_vals_sum,
                        curr_sample_sz, bin_step_size, max_extrap,
                        yield_vector);
      // no checking of curve in defect mode
      bootstrap_estimates.push_back(yield_vector);
      successful_bootstrap = true;
    }
    else {
      // refit curve for lower bound
      const ContinuedFractionApproximation lower_cfa(diagonal, max_terms);
      const ContinuedFraction lower_cf(
        lower_cfa.optimal_cont_frac_distinct(hist));

      // extrapolate the curve start
      if (lower_cf.is_valid()) {
        extrapolate_curve(lower_cf, initial_distinct, sample_vals_sum,
                          curr_sample_sz, bin_step_size, max_extrap,
                          yield_vector);
        // sanity check
        if (check_yield_estimates_stability(yield_vector)) {
          bootstrap_estimates.push_back(yield_vector);
          successful_bootstrap = true;
        }
      }
    }
    if (VERBOSE)
      cerr << (successful_bootstrap ? '.' : '_');
  }
  if (VERBOSE)
    cerr << endl;
  if (bootstrap_estimates.size() < n_bootstraps)
    throw runtime_error("too many defects in the approximation, "
                        "consider running in defect mode");
}

void
vector_median_and_ci(const vector<vector<double>> &bootstrap_estimates,
                     const double ci_level, vector<double> &yield_estimates,
                     vector<double> &lower_ci_lognorm,
                     vector<double> &upper_ci_lognorm) {
  yield_estimates.clear();
  lower_ci_lognorm.clear();
  upper_ci_lognorm.clear();
  assert(!bootstrap_estimates.empty());

  const size_t n_est = bootstrap_estimates.size();
  vector<double> estimates_row(n_est, 0.0);
  for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {
    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];

    double median_estimate, lower_ci_estimate, upper_ci_estimate;
    median_and_ci(estimates_row, ci_level, median_estimate, lower_ci_estimate,
                  upper_ci_estimate);
    std::sort(begin(estimates_row), end(estimates_row));

    yield_estimates.push_back(median_estimate);
    lower_ci_lognorm.push_back(lower_ci_estimate);
    upper_ci_lognorm.push_back(upper_ci_estimate);
  }
}

void
write_predicted_complexity_curve(const string &outfile, const double c_level,
                                 const double step_size,
                                 const vector<double> &yield_estimates,
                                 const vector<double> &yield_lower_ci_lognorm,
                                 const vector<double> &yield_upper_ci_lognorm) {
  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOWER_" << c_level << "CI\t"
      << "UPPER_" << c_level << "CI" << endl;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < yield_estimates.size(); ++i)
    out << (i + 1) * step_size << '\t' << yield_estimates[i] << '\t'
        << yield_lower_ci_lognorm[i] << '\t' << yield_upper_ci_lognorm[i]
        << endl;
}

// vals_hist[j] = n_{j} = # (counts = j)
// vals_hist_distinct_counts[k] = kth index j s.t. vals_hist[j] > 0
// stores kth index of vals_hist that is positive
// distinct_counts_hist[k] = vals_hist[vals_hist_distinct_counts[k]]
// stores the kth positive value of vals_hist
void
resample_hist(mt19937 &gen, const vector<size_t> &vals_hist_distinct_counts,
              const vector<double> &distinct_counts_hist,
              vector<double> &out_hist) {
  const size_t hist_size = distinct_counts_hist.size();
  vector<uint32_t> sample_distinct_counts_hist(hist_size, 0);

  const uint32_t distinct =
    accumulate(cbegin(distinct_counts_hist), cend(distinct_counts_hist), 0.0);

  multinomial(gen, distinct_counts_hist, distinct, sample_distinct_counts_hist);

  out_hist.clear();
  out_hist.resize(vals_hist_distinct_counts.back() + 1, 0.0);
  for (size_t i = 0; i < hist_size; i++)
    out_hist[vals_hist_distinct_counts[i]] = sample_distinct_counts_hist[i];
}

template <typename T>
T
median_from_sorted_vector(const vector<T> sorted_data, const size_t stride,
                          const size_t n) {
  if (n == 0 || sorted_data.empty())
    return 0.0;

  const size_t lhs = (n - 1) / 2;
  const size_t rhs = n / 2;

  if (lhs == rhs)
    return sorted_data[lhs * stride];

  return (sorted_data[lhs * stride] + sorted_data[rhs * stride]) / 2.0;
}

template <typename T>
T
quantile_from_sorted_vector(const vector<T> &sorted_data, const size_t stride,
                            const size_t n, const double f) {
  const double index = f * (n - 1);
  const size_t lhs = static_cast<int>(index);
  const double delta = index - lhs;

  if (n == 0 || sorted_data.empty())
    return 0.0;

  if (lhs == n - 1)
    return sorted_data[lhs * stride];

  return (1 - delta) * sorted_data[lhs * stride] +
         delta * sorted_data[(lhs + 1) * stride];
}

// Confidence interval stuff
void
median_and_ci(vector<double> estimates,  // by val so we can sort them
              const double ci_level, double &median_estimate,
              double &lower_ci_estimate, double &upper_ci_estimate) {
  assert(!estimates.empty());

  std::sort(begin(estimates), end(estimates));

  const double alpha = 1.0 - ci_level;
  const size_t N = estimates.size();

  median_estimate = median_from_sorted_vector(estimates, 1, N);
  lower_ci_estimate = quantile_from_sorted_vector(estimates, 1, N, alpha / 2);
  upper_ci_estimate =
    quantile_from_sorted_vector(estimates, 1, N, 1.0 - alpha / 2);
}
