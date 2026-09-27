// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#ifndef SRC_COMMON_HPP_
#define SRC_COMMON_HPP_

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <print>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

[[nodiscard]] auto
median_and_ci(const std::span<const double> values, const double ci_level)
  -> std::tuple<double, double, double>;

[[nodiscard]] auto
median_and_ci_md(const std::span<const std::vector<double>> values,
                 const double ci_level)
  -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

void
write_predicted_complexity_curve(
  const std::string &outfile,
  const double c_level,
  const double step_size,
  const std::vector<double> &yield_estimates,
  const std::vector<double> &yield_lower_ci_lognorm,
  const std::vector<double> &yield_upper_ci_lognorm);

template <typename T>
[[nodiscard]] auto
get_counts_from_hist(const std::vector<T> &h) -> T {
  T c = 0.0;
  for (auto i = 0u; i < std::size(h); ++i)
    c += i * h[i];
  return c;
}

template <typename uint_type>
void
multinomial(std::mt19937 &gen,
            const std::span<const double> mult_probs,
            uint_type trials,
            std::vector<uint_type> &result) {
  using binom_dist = std::binomial_distribution<uint_type>;

  result.clear();
  result.resize(std::size(mult_probs));

  double remaining_prob =
    std::reduce(std::cbegin(mult_probs), std::cend(mult_probs));

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

[[nodiscard]] auto
format_histogram(const auto &h) -> std::string {
  std::string s;
  for (auto i = 0u; i < std::size(h); ++i)
    if (h[i] > 0)
      s += std::format("{}\t{}", i, static_cast<std::uint32_t>(h[i]));
  return s;
}

auto
report_histogram(const std::string &outfile, const auto &h) -> void {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::print("{}", format_histogram(h));
}

#endif  // SRC_COMMON_HPP_
