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
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

[[nodiscard]] inline constexpr auto
positive_integer(const double x) -> bool {
  return std::nearbyint(x) == x;
}

[[nodiscard]] auto
median_and_ci(const std::span<const double> values, const double ci_level)
  -> std::tuple<double, double, double>;

[[nodiscard]] auto
median_and_ci_md(const std::span<const std::vector<double>> values,
                 const double ci_level)
  -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

void
write_complexity_curve(const std::string &outfile,
                       const std::span<const std::string> header,
                       const std::span<const double> points,
                       const std::span<const double> estimates,
                       const std::span<const double> lower_ci_lognorm,
                       const std::span<const double> upper_ci_lognorm);

void
write_complexity_curve(const std::string &outfile,
                       const std::span<const std::string> header,
                       const std::span<const double> points,
                       const std::span<const double> estimates);

template <typename T>
[[nodiscard]] auto
get_counts_from_hist(const std::vector<T> &h) -> T {
  T c = 0.0;
  for (auto i = 0U; i < std::size(h); ++i)
    c += i * h[i];
  return c;
}

[[nodiscard]] auto
format_histogram(const auto &h) -> std::string {
  std::string s;
  for (auto i = 0U; i < std::size(h); ++i)
    if (h[i] > 0)
      s += std::format("{}\t{}\n", i, static_cast<std::uint32_t>(h[i]));
  return s;
}

auto
report_histogram(const std::string &outfile, const auto &h) -> void {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::print(out, "{}", format_histogram(h));
}

#endif  // SRC_COMMON_HPP_
