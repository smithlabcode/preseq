// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "pop_size.hpp"

#include "cli_common.hpp"
#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include <CLI11/CLI11.hpp>
#include <libpreseq.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

auto
pop_size_main(const std::span<char *> args) -> int {
  static constexpr auto cmd_name = "pop_size";
  static constexpr auto max_iter_per_bootstrap = 100;
  static constexpr auto n_desired_steps_default = 50;
  static constexpr auto max_terms_default = 100;
  static constexpr auto default_max_extrap = 1'000'000'000LU;
  static constexpr auto min_required_counts = 4LU;
  static constexpr auto min_required_counts_error_message =
    "max count before zero is less than min required count (4)";

  const int argc = std::ssize(args);
  auto argv = std::data(args);

  std::string outfile;
  std::string infile;

  // NOLINTBEGIN(*-avoid-magic-numbers)
  std::uint32_t max_terms{max_terms_default};
  double max_extrap{};
  std::uint64_t n_desired_steps{n_desired_steps_default};
  std::uint64_t n_bootstraps = 100;
  constexpr int diagonal{};
  double c_level = 0.95;
  std::uint32_t seed = 408;
  // NOLINTEND(*-avoid-magic-numbers)

  // flags
  bool verbose{};
  bool single_estimate{};
  bool allow_defects{};

  std::uint32_t n_threads{1};

  CLI::App app{pop_size_about_msg};
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
    ->check(CLI::ExistingFile);
    // ->check(CLI::ReadPermission)
  app.add_option("-o,--output", outfile, "output file")
    ->option_text("FILE");
  app.add_option("-e,--extrap", max_extrap, "maximum extrapolation");
  app.add_option("-s,--steps", n_desired_steps, "number of steps");
  app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
  app.add_option("-c,--ci-level", c_level, "level for confidence intervals");
  app.add_option("-x,--terms", max_terms, "maximum terms in estimator");
  app.add_option("-r,--seed", seed, "seed for random number generator");
  app.add_flag("-q,--quick", single_estimate, "no bootstraps for confidence intervals");
  app.add_flag("-D,--defects", allow_defects, "no testing for defects");
  app.add_flag("-v,--verbose", verbose, "print more info");
  // clang-format on

  if (argc < 3) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  const auto input_format = get_input_format_type(infile);
  if (is_unknown(input_format)) {
    std::println("unknown input format");
    return EXIT_FAILURE;
  }

  if (verbose)
    std::println("INPUT FORMAT: {}", to_string(input_format));

  const auto [n_reads, counts_hist] = [&] {
    switch (input_format) {
    case input_format_type::hist:
      return load_histogram(infile);
    case input_format_type::counts:
      return load_counts(infile);
    case input_format_type::bam:
      return load_counts_BAM_se(n_threads, infile);
    default:  // case input_format_type::vals:
      return load_counts_bed_se(infile);
    }
  }();

  const std::uint64_t max_observed_count = std::size(counts_hist) - 1;
  const auto distinct_reads =
    std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

  // ensure that the max terms are acceptable
  const auto first_zero = preseq::get_first_zero(counts_hist);
  max_terms = std::min(max_terms, first_zero - 1);
  max_terms = max_terms - (max_terms % 2 == 1);
  if (max_extrap == 0.0)
    max_extrap = default_max_extrap * distinct_reads;

  const auto step_size =
    (max_extrap - distinct_reads) / static_cast<double>(n_desired_steps);

  const auto distinct_counts =
    std::ranges::count_if(counts_hist, [](const auto x) { return x > 0.0; });

  if (verbose)
    std::println("TOTAL READS     = {}"
                 "DISTINCT READS  = {}"
                 "DISTINCT COUNTS = {}"
                 "MAX COUNT       = {}"
                 "COUNTS OF 1     = {}"
                 "MAX TERMS       = {}",
                 n_reads,             //
                 distinct_reads,      //
                 distinct_counts,     //
                 distinct_counts,     //
                 max_observed_count,  //
                 counts_hist[1],      //
                 max_terms            //
    );

  // check to make sure library is not overly saturated
  const double two_fold_extrap = preseq::good_toulmin_2x(counts_hist);
  if (two_fold_extrap < 0.0)
    throw std::runtime_error(
      "Saturation expected at double initial sample size."
      " Unable to extrapolate");

  // check that min required count is satisfied
  if (max_terms < min_required_counts)
    throw std::runtime_error(min_required_counts_error_message);

  if (verbose)
    std::println("[ESTIMATING YIELD CURVE]");

  if (single_estimate) {
    const auto estimates = preseq::extrapolate_once(
      counts_hist, max_terms, diagonal, allow_defects, step_size, max_extrap);
    if (estimates.empty())
      throw std::runtime_error(
        "single estimate failed, run full mode for estimates");

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile.data());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    std::println(out, "TOTAL_READS\tEXPECTED_DISTINCT");
    std::println(out, "0\t0");
    for (auto i = 0LU; i < std::size(estimates); ++i)
      std::println(out, "{}\t{}",
                   std::round(static_cast<double>(i + 1) * step_size),
                   estimates[i]);
  }
  else {
    if (verbose)
      std::println("[BOOTSTRAPPING HISTOGRAM]");

    const std::uint64_t max_iter = max_iter_per_bootstrap * n_bootstraps;
    const auto cfa = preseq::cfa_t{
      .max_terms = max_terms,
      .diagonal = diagonal,
      .allow_defects = allow_defects,
      .rng_seed = seed,
      .target_bootstraps = n_bootstraps,
      .max_iterations = max_iter,
      .step_size = step_size,
      .max_extrap = max_extrap,
    };
    const auto estimates = extrapolate_bootstrap(cfa, counts_hist);
    if (verbose)
      std::println("[COMPUTING CONFIDENCE INTERVALS]");
    auto [medians, lower_ci_lognorms, upper_ci_lognorms] =
      median_and_ci_md(estimates, c_level);
    if (verbose)
      std::println("[WRITING OUTPUT]");

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    if (!outfile.empty() && !out)
      throw std::runtime_error("failed to open outfile file: " + outfile);

    const std::uint64_t n_ests = std::size(medians) - 1;
    if (n_ests < 2)
      throw std::runtime_error("problem with number of estimates in pop_size");

    const bool converged = (medians[n_ests] - medians[n_ests - 1] < 1.0);

    std::println(out, "pop_size_estimate\t"
                      "lower_ci\t"
                      "upper_ci");
    std::println(out, "{}\t{}\t{}", medians.back(), lower_ci_lognorms.back(),
                 upper_ci_lognorms.back());
    if (!converged)
      std::println(out, "\tnot_converged");
  }
  return EXIT_SUCCESS;
}
