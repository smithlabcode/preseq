// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "common.hpp"

#include <libpreseq.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

template <typename T>
[[nodiscard]] static auto
median_from_sorted(const std::span<T> &data, const auto n) -> T {
  static constexpr auto two = static_cast<T>(2);
  assert(std::ranges::is_sorted(data));
  if (n == 0 || data.empty())
    return 0.0;
  const auto lhs = (n - 1) / 2;
  const auto rhs = n / 2;
  return lhs == rhs ? data[lhs] : (data[lhs] + data[rhs]) / two;
}

template <typename T>
[[nodiscard]] static auto
quantile_from_sorted(const std::span<T> &data,
                     const std::uint64_t n,
                     const double f) -> T {
  assert(std::ranges::is_sorted(data));
  const auto index = f * static_cast<double>(n - 1);
  const std::size_t lhs = std::floor(index);
  const auto delta = index - lhs;
  if (n == 0 || data.empty())
    return 0.0;
  return static_cast<std::uint64_t>(lhs + 1) == n
           ? data[lhs]
           : (1 - delta) * data[lhs] + delta * data[lhs + 1];
}

/// Get median and confidence interval from sorted values
[[nodiscard]] auto
median_and_ci(const std::span<const double> values, const double ci_level)
  -> std::tuple<double, double, double> {
  assert(!values.empty() && std::ranges::is_sorted(values));
  const auto alpha = 1.0 - ci_level;
  const auto N = std::size(values);
  return std::tuple{median_from_sorted(values, N),
                    quantile_from_sorted(values, N, alpha / 2),
                    quantile_from_sorted(values, N, 1.0 - alpha / 2)};
}

/// Get median and confidence interval from each among multiple sets of values
[[nodiscard]] auto
median_and_ci_md(const std::span<const std::vector<double>> values,
                 const double ci_level)
  -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
  if (values.empty())
    throw std::runtime_error("empty input to median_and_ci");
  const auto N = std::size(values.front());
  const auto check_size = [N](const auto &x) { return std::size(x) == N; };
  if (!std::ranges::all_of(values, check_size))
    throw std::runtime_error("unequal vectors sizes in median_and_ci");

  std::vector<double> estimates;
  std::vector<double> lower;
  std::vector<double> upper;

  const auto n_est = std::ssize(values);
  for (auto i = 0U; i < N; ++i) {
    // Original data is not sorted; make a copy and sort
    auto estimates_row = std::vector(n_est, 0.0);
    for (auto k = 0; k < n_est; ++k)
      estimates_row[k] = values[k][i];
    std::ranges::sort(estimates_row);
    const auto est = median_and_ci(estimates_row, ci_level);
    estimates.push_back(std::get<0>(est));
    lower.push_back(std::get<1>(est));
    upper.push_back(std::get<2>(est));
  }
  return std::tuple{std::move(estimates), std::move(lower), std::move(upper)};
}

void
write_complexity_curve(const std::string &outfile,
                       const std::span<const std::string> header,
                       const std::span<const double> points,
                       const std::span<const double> estimates,
                       const std::span<const double> lower_ci,
                       const std::span<const double> upper_ci) {
  static constexpr auto row_fmt = "{}\t{}\t{}\t{}";
  if (std::size(header) != 4)
    throw std::runtime_error("incorrect header size (must be 4)");
  const auto n_points = std::size(points);
  if (n_points != std::size(estimates) || n_points != std::size(lower_ci) ||
      n_points != std::size(upper_ci))
    throw std::runtime_error(
      "unequal sized vectors in write_predicted_complexity_curve");
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::println(out, row_fmt, header[0], header[1], header[2], header[3]);
  const auto toi = [](const auto x) -> std::int64_t { return std::round(x); };
  for (auto i = 0; i < std::ssize(estimates); ++i)
    std::println(out, row_fmt,       //
                 toi(points[i]),     //
                 toi(estimates[i]),  //
                 toi(lower_ci[i]),   //
                 toi(upper_ci[i]));
}

void
write_complexity_curve(const std::string &outfile,
                       const std::span<const std::string> header,
                       const std::span<const double> points,
                       const std::span<const double> estimates) {
  static constexpr auto row_fmt = "{}\t{}";
  if (std::size(header) != 2)
    throw std::runtime_error("incorrect header size (must be 2)");
  if (std::size(points) != std::size(estimates))
    throw std::runtime_error("unequal sized vectors in write_complexity_curve");
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::println(out, row_fmt, header[0], header[1]);
  const auto toi = [](const auto x) -> std::int64_t { return std::round(x); };
  for (auto i = 0ul; i < std::size(estimates); ++i)
    std::println(out, row_fmt, toi(points[i]), toi(estimates[i]));
}
