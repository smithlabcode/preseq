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

#ifndef SRC_COMMON_HPP_
#define SRC_COMMON_HPP_

#include <cstddef>  // std::size_t
#include <cstdint>  // std::uint64_t
#include <random>
#include <stdexcept>
#include <string>
#include <vector>  // has std::size

double
GoodToulmin2xExtrap(const std::vector<double> &counts_hist);

double
interpolate_distinct(const std::vector<double> &hist, const std::size_t N,
                     const std::size_t S, const std::size_t n);

bool
extrap_single_estimate(const bool VERBOSE, const bool allow_defects,
                       const std::vector<double> &hist, std::size_t max_terms,
                       const int diagonal, const double step_size,
                       const double max_extrap,
                       std::vector<double> &yield_estimate);

void
extrap_bootstrap(const bool VERBOSE, const bool allow_defects,
                 const std::uint64_t seed, const std::vector<double> &orig_hist,
                 const std::size_t n_bootstraps,
                 const std::size_t orig_max_terms, const int diagonal,
                 const double bin_step_size, const double max_extrap,
                 const std::size_t max_iter,
                 std::vector<std::vector<double>> &bootstrap_estimates);

void
vector_median_and_ci(
  const std::vector<std::vector<double>> &bootstrap_estimates,
  const double ci_level, std::vector<double> &yield_estimates,
  std::vector<double> &lower_ci_lognorm, std::vector<double> &upper_ci_lognorm);

void
write_predicted_complexity_curve(
  const std::string &outfile, const double c_level, const double step_size,
  const std::vector<double> &yield_estimates,
  const std::vector<double> &yield_lower_ci_lognorm,
  const std::vector<double> &yield_upper_ci_lognorm);

template <typename T>
T
get_counts_from_hist(const std::vector<T> &h) {
  T c = 0.0;
  for (auto i = 0u; i < std::size(h); ++i)
    c += i * h[i];
  return c;
}

double
factorial(double x);

void
resample_hist(std::mt19937 &gen,
              const std::vector<std::size_t> &vals_hist_distinct_counts,
              const std::vector<double> &distinct_counts_hist,
              std::vector<double> &out_hist);

void
median_and_ci(std::vector<double> estimates,  // by val so we can sort them
              const double ci_level, double &median_estimate,
              double &lower_ci_estimate, double &upper_ci_estimate);

template <typename uint_type>
void
multinomial(std::mt19937 &gen, const std::vector<double> &mult_probs,
            uint_type trials, std::vector<uint_type> &result) {
  typedef std::binomial_distribution<uint32_t> binom_dist;

  result.clear();
  result.resize(std::size(mult_probs));

  double remaining_prob =
    std::accumulate(std::begin(mult_probs), std::end(mult_probs), 0.0);

  auto r = std::begin(result);
  auto p = std::begin(mult_probs);

  while (p != std::end(mult_probs)) {  // iterate to sample for each category
    *r = binom_dist(trials, (*p) / remaining_prob)(gen);  // take the sample

    remaining_prob -= *p++;  // update remaining probability mass
    trials -= *r++;          // update remaining trials needed
  }

  if (trials > 0)
    throw std::runtime_error("multinomial sampling failed");
}

#endif  // SRC_COMMON_HPP_
