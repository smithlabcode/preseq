// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "bound_pop.hpp"

#include "cli_common.hpp"
#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include <libpreseq.hpp>

#include <CLI11/CLI11.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>  // IWYU pragma: keep
#include <numeric>
#include <print>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

// bounding n_0
auto
bound_pop_main(const std::span<char *> args) -> int {
  static constexpr auto cmd_name = "bound_pop";
  const auto normalize =
    [](auto &x) {  // cppcheck-suppress[constParameterReference]
      const auto d = std::reduce(std::cbegin(x), std::cend(x));
      std::ranges::transform(x, std::begin(x),
                             [&](const auto y) { return y / d; });
    };

  const int argc = std::ssize(args);
  auto argv = std::data(args);

  bool verbose{};
  bool quick_mode{};

  std::string infile;
  std::string outfile;

  // NOLINTBEGIN(*-avoid-magic-numbers)
  std::size_t max_num_points = 10;
  double tolerance = 1e-20;
  std::size_t n_bootstraps = 500;
  double c_level = 0.95;
  constexpr std::size_t max_iter = 100;
  std::uint32_t seed = 408;
  // NOLINTEND(*-avoid-magic-numbers)

  std::uint32_t n_threads{1};

  CLI::App app{bound_pop_about_msg};
  argv = app.ensure_utf8(argv);
  app.failure_message(
    [&](const CLI::App *a, const CLI::Error &b) -> std::string {
      return std::format("preseq command: {}\n", cmd_name) +
             CLI::FailureMessage::simple(a, b);
    });
  app.usage(std::format("Usage: preseq {} [OPTIONS]", cmd_name));
  if (argc >= 2)
    app.footer(" ");  // DESCRIPTION

  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  app.formatter(std::make_shared<preseq_formatter>());
  app.get_formatter()->column_width(24);
  app.get_formatter()->long_option_alignment_ratio(0.2);
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  // clang-format off
  app.add_option("INPUT", infile, "input file (hist/BED/BAM/SAM/values)")
    ->required()
    ->option_text(" ")
    ->check(CLI::ExistingFile)
    ->check(CLI::ReadPermissions);
  app.add_option("-o,--output", outfile, "output file")
    ->option_text("FILE")
    ->check(CLI::WritePermissions);
  app.add_option("-m,--max-points", max_num_points,
                 "max points in estimates by quadrature")
    ->option_text("INT")
    ->check(CLI::PositiveNumber);
  // app.add_option("-t,--tol", tolerance, "numerical tolerance");
  app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
  app.add_option("-c,--cval", c_level, "level for confidence intervals")
    ->ignore_underscore();
  app.add_option("-r,--seed", seed, "seed for random number generator");
  app.add_flag("-q,--quick", quick_mode, "no bootstraps when making estimates");
  app.add_flag("-v,--verbose", verbose, "print moments and boostraps with output");
  // clang-format on

  if (argc < 2) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  const auto input_format = get_input_format_type(infile);
  if (is_unknown(input_format)) {
    std::println("unknown input format");
    return EXIT_FAILURE;
  }

  const auto [n_obs, counts_hist] = [&] {
    if (is_hist(input_format))
      return load_histogram(infile);
    if (is_counts(input_format))
      return load_counts(infile);
    if (is_bam(input_format))
      return load_counts_BAM_se(n_threads, infile);
    //  if (is_bed(input_format))
    return load_counts_bed_se(infile);
  }();

  const double distinct_obs =
    std::reduce(std::cbegin(counts_hist), std::cend(counts_hist));

  auto measure_moments = preseq::get_measure_moments(counts_hist);

  std::vector<nlohmann::json> bootstraps;
  nlohmann::json output;

  if (quick_mode) {
    if (std::size(measure_moments) > 2 * max_num_points)
      measure_moments.resize(2 * max_num_points);

    auto [n_points, pd_mom_seq] =
      preseq::get_posdef_moment_sequence(measure_moments, tolerance);
    preseq::moment_sequence_t obs_mom_seq(pd_mom_seq);

    auto [points, weights] =
      obs_mom_seq.lower_quadrature_rules(n_points, tolerance, max_iter);

    normalize(weights);

    const auto n_1 = counts_hist[1];
    const auto term = [n_1](const auto w, const auto p) { return n_1 * w / p; };
    auto estimated_unobs =
      std::inner_product(std::cbegin(weights), std::cend(weights),
                         std::cbegin(points), 0.0, std::plus<>(), term);
    estimated_unobs = std::max(estimated_unobs, 0.0) + distinct_obs;
    if (estimated_unobs == distinct_obs)
      n_points = 0;

    output = {
      {"quadrature_estimated_unobs", estimated_unobs},
      {"n_points", n_points},
    };
  }
  else {
    // do bootstraps
    std::mt19937 gen{seed};
    std::vector<double> quad_estimates;
    const auto lim = std::min(max_iter, n_bootstraps);
    for (auto i = 0u; i < lim; ++i) {
      auto sample_hist = preseq::resample_hist(gen, counts_hist);
      const double sampled_distinct =
        std::reduce(std::cbegin(sample_hist), std::cend(sample_hist));

      // initialize moments, 0-th moment is 1
      std::vector<double> bootstrap_moments(1, 1.0);
      // moments[r] = (r + 1)! n_{r+1} / n_1
      for (std::size_t j = 0; j < 2 * max_num_points; ++j)
        bootstrap_moments.push_back(std::exp(lnfact(j + 3) +
                                             std::log(sample_hist[j + 2]) -
                                             std::log(sample_hist[1])));
      auto [n_points, pd_mom_seq] =
        preseq::get_posdef_moment_sequence(bootstrap_moments, tolerance);

      preseq::moment_sequence_t obs_mom_seq(pd_mom_seq);
      auto [points, weights] =
        obs_mom_seq.lower_quadrature_rules(n_points, tolerance, max_iter);
      normalize(weights);

      const auto n_1 = counts_hist[1];
      const auto term = [n_1](const auto w, const auto p) {
        return n_1 * w / p;
      };
      auto estimated_unobs =
        std::inner_product(std::cbegin(weights), std::cend(weights),
                           std::cbegin(points), 0.0, std::plus{}, term);
      estimated_unobs = std::max(estimated_unobs, 0.0) + sampled_distinct;
      if (verbose)
        bootstraps.push_back(nlohmann::json({
          {"bootstrapped_moments", bootstrap_moments},
          {"moment_sequence", to_string(obs_mom_seq)},
          {"points", points},
          {"weights", weights},
          {"estimated_unobs", estimated_unobs},
        }));
      quad_estimates.push_back(estimated_unobs);
    }
    const auto [median, lower_ci, upper_ci] =
      median_and_ci(quad_estimates, c_level);
    output = {
      {"median_estimated_unobs", median},
      {"lower_ci", lower_ci},
      {"upper_ci", upper_ci},
    };
  }

  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!outfile.empty() && !out)
    throw std::runtime_error("failed to open output file: " + outfile);

  output["total_observations"] = n_obs;
  output["distinct_observations"] = distinct_obs;
  output["max_count"] = std::size(counts_hist) - 1;
  output["observed_moments"] = measure_moments;
  if (verbose && !quick_mode)
    output["bootstraps"] = bootstraps;
  std::println(out, "{}", output.dump(4));
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
