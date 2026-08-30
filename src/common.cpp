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
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "common.hpp"

#include "continued_fraction.hpp"
#include "lnfact.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

[[nodiscard]] auto
GoodToulmin2xExtrap(const std::vector<double> &counts_hist) -> double {
  double two_fold_extrap = 0.0;
  for (std::size_t i = 0; i < std::size(counts_hist); ++i) {
    const int sign = (i % 2 == 0) ? -1 : 1;  // (-1)^(i+1)
    two_fold_extrap += sign * counts_hist[i];
  }
  return two_fold_extrap;
}

// interpolate by explicit calculating the expectation for sampling without
// replacement; see K.L Heck 1975
//
// -- N total sample size; S the total number of distincts
// -- n sub sample size
[[nodiscard]] auto
interpolate_distinct(const std::vector<double> &hist, const std::size_t N,
                     const std::size_t S, const std::size_t n) -> double {
  const double log_denom = lnfact(N + 1) - lnfact(n + 1) - lnfact(N - n + 1);
  std::vector<double> numerator(std::size(hist), 0);
  for (std::size_t i = 1; i < std::size(hist); ++i) {
    // N - i - n + 1 should be greater than 0
    if (N < i + n)
      continue;
    const auto x = lnfact(N - i + 1) - lnfact(n + 1) - lnfact(N - i - n + 1);
    numerator[i] = std::exp(x - log_denom) * hist[i];
  }
  return S - std::reduce(std::cbegin(numerator), std::cend(numerator));
}

static auto
extrapolate_curve(const ContinuedFraction &the_cf,
                  const double initial_distinct, const double vals_sum,
                  const double initial_sample_size, const double step_size,
                  const double max_sample_size,
                  std::vector<double> &estimates) {
  double curr_samp_sz = initial_sample_size;
  while (curr_samp_sz < max_sample_size) {
    const double fold = (curr_samp_sz - vals_sum) / vals_sum;
    assert(fold >= 0.0);
    estimates.push_back(initial_distinct + fold * the_cf(fold));
    curr_samp_sz += step_size;
  }
}

[[nodiscard]] auto
extrap_single_estimate(const bool VERBOSE, const bool allow_defects,
                       const std::vector<double> &hist, std::size_t max_terms,
                       const int diagonal, const double step_size,
                       const double max_extrap,
                       std::vector<double> &yield_estimate) -> bool {
  yield_estimate.clear();

  const double vals_sum = get_counts_from_hist(hist);
  const double initial_distinct =
    std::accumulate(std::cbegin(hist), std::cend(hist), 0.0);

  // interpolate complexity curve by random sampling w/out replacement
  const std::size_t upper_limit = vals_sum;
  const std::size_t step = step_size;
  auto sample = static_cast<std::size_t>(step_size);
  for (; sample < upper_limit; sample += step)
    yield_estimate.push_back(
      interpolate_distinct(hist, upper_limit, initial_distinct, sample));

  // ensure that the max terms are acceptable
  std::size_t first_zero = 1;
  while (first_zero < std::size(hist) && hist[first_zero] > 0)
    ++first_zero;

  // Ensure we are not using a zero term
  max_terms = std::min(max_terms, first_zero - 1);

  // refit curve for lower bound (degree of approx is 1 less than
  // max_terms)
  max_terms = max_terms - (max_terms % 2 == 1);

  if (allow_defects) {
    std::vector<double> ps_coeffs;
    for (std::size_t j = 1; j <= max_terms; ++j)
      ps_coeffs.push_back(hist[j] * std::pow(-1.0, j + 1));

    const ContinuedFraction defect_cf(ps_coeffs, diagonal, max_terms);

    extrapolate_curve(defect_cf, initial_distinct, vals_sum, sample, step_size,
                      max_extrap, yield_estimate);

    if (VERBOSE)
      std::println(std::cerr, "{}", defect_cf);
    // NO FAIL! defect mode doesn't care about failure
  }
  else {
    const ContinuedFractionApproximation lower_cfa(diagonal, max_terms);
    const ContinuedFraction lower_cf(lower_cfa.optimal_cf_distinct(hist));

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
      std::cerr << lower_cf << '\n';
  }
  // SUCCESS!!
  return true;
}

void
extrap_bootstrap(const bool VERBOSE, const bool allow_defects,
                 const std::uint32_t rng_seed,
                 const std::vector<double> &orig_hist,
                 const std::size_t n_bootstraps,
                 const std::size_t orig_max_terms, const int diagonal,
                 const double bin_step_size, const double max_extrap,
                 const std::size_t max_iter,
                 std::vector<std::vector<double>> &bootstrap_estimates) {
  static constexpr auto progress_width = 72;  // bootstrap success progress
  // clear returning vectors
  bootstrap_estimates.clear();

  // setup rng
  std::mt19937 rng(rng_seed);

  const double initial_distinct =
    std::accumulate(std::cbegin(orig_hist), std::cend(orig_hist), 0.0);

  std::vector<std::size_t> orig_hist_distinct_counts;
  std::vector<double> distinct_orig_hist;
  for (std::size_t i = 0; i < std::size(orig_hist); ++i)
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }

  for (std::size_t iter = 0;
       (iter < max_iter && std::size(bootstrap_estimates) < n_bootstraps);
       ++iter) {
    if (VERBOSE && iter > 0 && iter % progress_width == 0)
      std::cerr << '\n';

    std::vector<double> yield_vector;
    std::vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);

    const double sample_vals_sum = get_counts_from_hist(hist);

    // resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    // compute complexity curve by random sampling w/out replacement
    const std::size_t distinct =
      std::accumulate(std::cbegin(hist), std::cend(hist), 0.0);
    std::size_t curr_sample_sz = bin_step_size;
    while (curr_sample_sz < sample_vals_sum) {
      yield_vector.push_back(
        interpolate_distinct(hist, sample_vals_sum, distinct, curr_sample_sz));
      curr_sample_sz += bin_step_size;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::size_t first_zero = 1;
    while (first_zero < std::size(hist) && hist[first_zero] > 0)
      ++first_zero;

    std::size_t max_terms = std::min(orig_max_terms, first_zero - 1);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 1);

    bool successful_bootstrap = false;
    // defect mode, simple extrapolation
    if (allow_defects) {
      std::vector<double> ps_coeffs;
      for (std::size_t j = 1; j <= max_terms; ++j)
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
      const ContinuedFraction lower_cf(lower_cfa.optimal_cf_distinct(hist));

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
      std::cerr << (successful_bootstrap ? '.' : '_');
  }
  if (VERBOSE)
    std::cerr << '\n';
  if (std::size(bootstrap_estimates) < n_bootstraps)
    throw std::runtime_error("too many defects in the approximation, "
                             "consider running in defect mode");
}

void
vector_median_and_ci(
  const std::vector<std::vector<double>> &bootstrap_estimates,
  const double ci_level, std::vector<double> &yield_estimates,
  std::vector<double> &lower_ci_lognorm,
  std::vector<double> &upper_ci_lognorm) {
  yield_estimates.clear();
  lower_ci_lognorm.clear();
  upper_ci_lognorm.clear();
  assert(!bootstrap_estimates.empty());

  const std::size_t n_est = std::size(bootstrap_estimates);
  std::vector<double> estimates_row(n_est, 0.0);
  for (std::size_t i = 0; i < std::size(bootstrap_estimates[0]); ++i) {
    // estimates is in wrong order, work locally on const val
    for (std::size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];

    double median_estimate{};
    double lower_ci_estimate{};
    double upper_ci_estimate{};
    median_and_ci(estimates_row, ci_level, median_estimate, lower_ci_estimate,
                  upper_ci_estimate);
    std::sort(std::begin(estimates_row), std::end(estimates_row));

    yield_estimates.push_back(median_estimate);
    lower_ci_lognorm.push_back(lower_ci_estimate);
    upper_ci_lognorm.push_back(upper_ci_estimate);
  }
}

void
write_predicted_complexity_curve(
  const std::string &outfile, const double c_level, const double step_size,
  const std::vector<double> &yield_estimates,
  const std::vector<double> &yield_lower_ci_lognorm,
  const std::vector<double> &yield_upper_ci_lognorm) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);

  // clang-format off
  std::println(out, "TOTAL_READS\t"
               "EXPECTED_DISTINCT\t"
               "LOWER_{0}CI\t"
               "UPPER_{0}CI",
               c_level);
  // clang-format on

  std::println(out, "0\t0\t0\t0");
  for (auto i = 0ul; i < std::size(yield_estimates); ++i)
    std::println(out, "{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}",  //
                 (i + 1) * step_size,                    //
                 yield_estimates[i],                     //
                 yield_lower_ci_lognorm[i],              //
                 yield_upper_ci_lognorm[i]               //
    );
}

// vals_hist[j] = n_{j} = # (counts = j)
// vals_hist_distinct_counts[k] = kth index j s.t. vals_hist[j] > 0
// stores kth index of vals_hist that is positive
// distinct_counts_hist[k] = vals_hist[vals_hist_distinct_counts[k]]
// stores the kth positive value of vals_hist
void
resample_hist(std::mt19937 &gen,
              const std::vector<std::size_t> &vals_hist_distinct_counts,
              const std::vector<double> &distinct_counts_hist,
              std::vector<double> &out_hist) {
  assert(std::ranges::is_sorted(vals_hist_distinct_counts));
  const std::size_t hist_size = std::size(distinct_counts_hist);
  std::vector<std::uint64_t> sample_distinct_counts_hist(hist_size, 0);
  const std::uint64_t distinct = std::reduce(std::cbegin(distinct_counts_hist),
                                             std::cend(distinct_counts_hist));
  multinomial(gen, distinct_counts_hist, distinct, sample_distinct_counts_hist);
  out_hist.resize(vals_hist_distinct_counts.back() + 1, 0.0);
  for (std::size_t i = 0; i < hist_size; ++i)
    out_hist[vals_hist_distinct_counts[i]] = sample_distinct_counts_hist[i];
}

template <typename T>
[[nodiscard]] auto
median_from_sorted_vector(const std::vector<T> &sorted_data,
                          const std::size_t n) -> T {
  if (n == 0 || sorted_data.empty())
    return 0.0;

  const std::size_t lhs = (n - 1) / 2;
  const std::size_t rhs = n / 2;

  if (lhs == rhs)
    return sorted_data[lhs];

  return (sorted_data[lhs] + sorted_data[rhs]) / static_cast<T>(2);
}

template <typename T>
[[nodiscard]] auto
quantile_from_sorted_vector(const std::vector<T> &sorted_data,
                            const std::size_t n, const double f) -> T {
  const double index = f * (n - 1);
  const std::size_t lhs = static_cast<std::size_t>(index);
  const double delta = index - lhs;

  if (n == 0 || sorted_data.empty())
    return 0.0;

  if (lhs + 1 == n)
    return sorted_data[lhs];

  return (1 - delta) * sorted_data[lhs] + delta * sorted_data[lhs + 1];
}

// Confidence interval stuff
void
median_and_ci(std::vector<double> estimates,  // by val so we can sort them
              const double ci_level, double &median_estimate,
              double &lower_ci_estimate, double &upper_ci_estimate) {
  assert(!estimates.empty());
  std::ranges::sort(estimates);

  const double alpha = 1.0 - ci_level;
  const std::size_t N = std::size(estimates);

  median_estimate = median_from_sorted_vector(estimates, N);
  lower_ci_estimate = quantile_from_sorted_vector(estimates, N, alpha / 2);
  upper_ci_estimate =
    quantile_from_sorted_vector(estimates, N, 1.0 - alpha / 2);
}

// NOLINTEND(*-narrowing-conversions)
