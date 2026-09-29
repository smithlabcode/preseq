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
median_from_sorted(const std::span<T> &sorted_data, const auto n) -> T {
  assert(std::ranges::is_sorted(sorted_data));
  if (n == 0 || sorted_data.empty())
    return 0.0;

  const auto lhs = (n - 1) / 2;
  const auto rhs = n / 2;

  if (lhs == rhs)
    return sorted_data[lhs];

  return (sorted_data[lhs] + sorted_data[rhs]) / static_cast<T>(2);
}

template <typename T>
[[nodiscard]] static auto
quantile_from_sorted(const std::span<T> &sorted_data,
                     const std::uint64_t n,
                     const double f) -> T {
  assert(std::ranges::is_sorted(sorted_data));
  const auto index = f * static_cast<double>(n - 1);
  const auto lhs = std::floor(index);
  const auto delta = index - lhs;
  if (n == 0 || sorted_data.empty())
    return 0.0;
  if (static_cast<std::uint64_t>(lhs + 1) == n)
    return sorted_data[lhs];
  return (1 - delta) * sorted_data[lhs] + delta * sorted_data[lhs + 1];
}

// Confidence interval stuff
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

[[nodiscard]] auto
median_and_ci_md(const std::span<const std::vector<double>> values,
                 const double ci_level)
  -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
  if (values.empty())
    throw std::runtime_error("empty input to median_and_ci");
  const auto N = std::size(values.front());
  if (!std::ranges::all_of(values,
                           [N](const auto &x) { return std::size(x) == N; }))
    throw std::runtime_error("unequal vectors sizes in median_and_ci");

  std::vector<double> estimates;
  std::vector<double> lower_ci;
  std::vector<double> upper_ci;

  const auto n_est = std::ssize(values);
  for (auto i = 0U; i < N; ++i) {
    // estimates is in wrong order, work locally on const val
    auto estimates_row = std::vector(n_est, 0.0);
    for (auto k = 0; k < n_est; ++k)
      estimates_row[k] = values[k][i];

    std::ranges::sort(estimates_row);
    const auto est = median_and_ci(estimates_row, ci_level);

    estimates.push_back(std::get<0>(est));
    lower_ci.push_back(std::get<1>(est));
    upper_ci.push_back(std::get<2>(est));
  }
  return std::tuple{std::move(estimates), std::move(lower_ci),
                    std::move(upper_ci)};
}

void
write_complexity_curve(const std::string &outfile,
                       const std::span<const std::string> header,
                       const std::span<const double> points,
                       const std::span<const double> estimates,
                       const std::span<const double> lower_ci_lognorm,
                       const std::span<const double> upper_ci_lognorm) {
  if (std::size(header) != 4)
    throw std::runtime_error("incorrect header size (must be 4)");
  if (std::size(points) != std::size(estimates) ||
      std::size(points) != std::size(lower_ci_lognorm) ||
      std::size(points) != std::size(upper_ci_lognorm))
    throw std::runtime_error(
      "unequal sized vectors in write_predicted_complexity_curve");
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::println(out, "{}\t{}\t{}\t{}", header[0], header[1], header[2],
               header[3]);
  for (auto i = 0ul; i < std::size(estimates); ++i)
    std::println(
      out, "{}\t{}\t{}\t{}",                                        //
      static_cast<std::uint64_t>(std::round(points[i])),            //
      static_cast<std::uint64_t>(std::round(estimates[i])),         //
      static_cast<std::uint64_t>(std::round(lower_ci_lognorm[i])),  //
      static_cast<std::uint64_t>(std::round(upper_ci_lognorm[i]))   //
    );
}

void
write_complexity_curve(const std::string &outfile,
                       const std::span<const std::string> header,
                       const std::span<const double> points,
                       const std::span<const double> estimates) {
  if (std::size(header) != 2)
    throw std::runtime_error("incorrect header size (must be 2)");
  if (std::size(points) != std::size(estimates))
    throw std::runtime_error(
      "unequal sized vectors in write_predicted_complexity_curve");
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::println(out, "{}\t{}", header[0], header[1]);
  for (auto i = 0ul; i < std::size(estimates); ++i)
    std::println(out, "{}\t{}",                                      //
                 static_cast<std::uint64_t>(std::round(points[i])),  //
                 static_cast<std::uint64_t>(std::round(estimates[i])));
}
