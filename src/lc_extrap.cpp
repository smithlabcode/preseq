// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "lc_extrap.hpp"

#include "cli_common.hpp"
#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include <CLI11/CLI11.hpp>
#include <libpreseq.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <numeric>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

[[nodiscard]] static auto
get_points(const auto n, const auto s) {
  return std::views::iota(0, n) |
         std::views::transform([&](const auto x) { return (x + 1) * s; }) |
         std::ranges::to<std::vector>();
}

auto
lc_extrap_main(const std::span<char *> args) -> int {
  using std::string_literals::operator""s;
  static constexpr auto cmd_name = "lc_extrap";
  static constexpr auto max_extrap_default = 10'000'000'000.0;
  static const std::size_t min_required_counts = 4;
  static const std::string min_required_counts_error_message =
    std::format("max count before zero is less than min required count ({}) "
                "duplicates removed",
                min_required_counts);

  const int argc = std::ssize(args);
  auto argv = std::data(args);

  std::string outfile;
  std::string input_file_name;

  // NOLINTBEGIN(*-avoid-magic-numbers)
  constexpr int diagonal{0};
  std::uint32_t orig_max_terms{100};
  double max_extrap{max_extrap_default};
  double step_size{1e6};
  std::size_t n_bootstraps{100};
  double c_level{0.95};
  std::uint32_t seed{408};
  // NOLINTEND(*-avoid-magic-numbers)

  bool verbose{};
  bool single_estimate{};
  bool allow_defects{};

  std::uint32_t n_threads{1};

  CLI::App app{lc_extrap_about_msg};
  argv = app.ensure_utf8(argv);
  app.failure_message(
    [&](const CLI::App *a, const CLI::Error &b) -> std::string {
      return std::format("preseq command: {}\n", cmd_name) +
             CLI::FailureMessage::simple(a, b);
    });
  app.usage(std::format("Usage: preseq {} [OPTIONS]", cmd_name));
  if (argc >= 2)
    app.footer(lc_extrap_footer_msg);

  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  app.formatter(std::make_shared<preseq_formatter>());
  app.get_formatter()->column_width(22);
  app.get_formatter()->long_option_alignment_ratio(0.2);
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  // clang-format off
  app.require_subcommand(0);
  app.add_option("INPUT", input_file_name, "input file (hist/BED/BAM/SAM/values)")
    ->required()
    ->option_text(" ")
    ->check(CLI::ExistingFile);
  // ->check(CLI::ReadPermission);
  app.add_option("-o,--output", outfile, "output filename (directory must exist)")
    ->option_text("FILE")
    ->required();
  // ->check(CLI::WritePermission)
  app.add_option("-e,--extrap", max_extrap, "maximum extrapolation (must be integer)")
    ->option_text("INT")
    ->check(CLI::PositiveNumber);
  app.add_option("-s,--step", step_size, "extrapolation step size");
  app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
  app.add_option("-c,--cval", c_level, "level for confidence intervals");
  app.add_option("-x,--terms", orig_max_terms, "maximum terms in estimator");
  app.add_option("-r,--seed", seed, "seed for random number generator");
  app.add_flag("-q,--quick", single_estimate, "no bootstraps for confidence intervals");
  app.add_flag("-D,--defects", allow_defects, "no testing for defects");
  app.add_flag("-v,--verbose", verbose, "print more info");
  app.set_help_flag("-h,--help", "print a detailed help message and exit");
  // clang-format on

  if (argc < 2) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  const auto input_format = get_input_format_type(input_file_name);
  if (is_unknown(input_format)) {
    std::println("unknown input format");
    return EXIT_FAILURE;
  }

  if (verbose)
    std::println("INPUT FORMAT: {}", to_string(input_format));

  const auto [n_reads, counts_hist] = [&] {
    switch (input_format) {
    case input_format_type::hist:
      return load_histogram(input_file_name);
    case input_format_type::counts:
      return load_counts(input_file_name);
    case input_format_type::bam:
      return load_counts_BAM_se(n_threads, input_file_name);
    default:  // case input_format_type::vals:
      return load_counts_bed_se(input_file_name);
    }
  }();

  const auto max_observed_count = std::size(counts_hist) - 1;
  const double distinct_reads =
    std::reduce(std::cbegin(counts_hist), std::cend(counts_hist));

  // make sure the max terms is at most one less than the first zero
  const auto first_zero = preseq::get_first_zero(counts_hist);
  orig_max_terms = std::min(orig_max_terms, first_zero - 1);
  orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

  const auto distinct_counts =
    std::ranges::count_if(counts_hist, [](const double x) { return x > 0.0; });

  if (verbose)
    std::print("TOTAL READS     = {}\n"
               "DISTINCT READS  = {}\n"
               "DISTINCT COUNTS = {}\n"
               "MAX COUNT       = {}\n"
               "COUNTS OF 1     = {}\n"
               "MAX TERMS       = {}\n",
               n_reads,             //
               distinct_reads,      //
               distinct_counts,     //
               max_observed_count,  //
               counts_hist[1],      //
               orig_max_terms);

  // check to make sure library is not overly saturated
  const double two_fold_extrap = preseq::good_toulmin_2x(counts_hist);
  if (two_fold_extrap < 0.0)
    throw std::runtime_error(
      "Saturation expected at double initial sample size. "
      "Unable to extrapolate.");

  // check that min required count is satisfied
  if (orig_max_terms < min_required_counts)
    throw std::runtime_error(min_required_counts_error_message);

  if (verbose)
    std::println("[ESTIMATING YIELD CURVE]");

  if (single_estimate) {
    auto estimates =
      preseq::extrapolate_once(counts_hist, orig_max_terms, diagonal,
                               allow_defects, step_size, max_extrap);
    if (estimates.empty())
      throw std::runtime_error("single estimate failed; run in full mode");

    auto points = get_points(std::ssize(estimates), step_size);

    points.insert(std::cbegin(points), 0.0);
    estimates.insert(std::cbegin(estimates), 0.0);

    const auto header = std::vector({
      "TOTAL_READS"s,
      "EXPECTED_DISTINCT"s,
    });

    write_complexity_curve(outfile, header, points, estimates);
  }
  else {
    if (verbose)
      std::println("[BOOTSTRAPPING HISTOGRAM]");
    const std::size_t max_iter = 100 * n_bootstraps;
    const auto cfa = preseq::cfa_t{
      .max_terms = orig_max_terms,
      .diagonal = diagonal,
      .allow_defects = allow_defects,
      .rng_seed = seed,
      .target_bootstraps = n_bootstraps,
      .max_iterations = max_iter,
      .step_size = step_size,
      .max_extrap = max_extrap,
    };
    auto bootstrap_estimates = preseq::extrapolate_bootstrap(cfa, counts_hist);
    if (verbose)
      std::println("[COMPUTING CONFIDENCE INTERVALS]");
    auto [medians, lower_ci, upper_ci] =
      median_and_ci_md(bootstrap_estimates, c_level);
    if (verbose)
      std::println("[WRITING OUTPUT]");

    auto points = get_points(std::ssize(medians), step_size);

    points.insert(std::cbegin(points), 0.0);
    medians.insert(std::cbegin(medians), 0.0);
    lower_ci.insert(std::cbegin(lower_ci), 0.0);
    upper_ci.insert(std::cbegin(upper_ci), 0.0);

    const auto header = std::vector({
      "TOTAL_READS"s,
      "EXPECTED_DISTINCT"s,
      std::format("LOWER_{}CI", c_level),
      std::format("UPPER_{}CI", c_level),
    });

    write_complexity_curve(outfile, header, points, medians, lower_ci,
                           upper_ci);
  }

  return EXIT_SUCCESS;
}
