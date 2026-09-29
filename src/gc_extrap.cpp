// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "gc_extrap.hpp"

#include "cli_common.hpp"
#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include <CLI11/CLI11.hpp>
#include <libpreseq.hpp>

#include <algorithm>
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

static auto
write_output(const std::string &outfile,          //
             const std::uint32_t bin_size,        //
             const std::uint32_t base_step_size,  //
             const std::vector<double> &coverage_estimates) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::println(out, "TOTAL_BASES\tEXPECTED_DISTINCT");
  std::println(out, "0\t0");
  for (auto i = 0LU; i < std::size(coverage_estimates); ++i)
    std::println(out, "{}\t{}",
                 static_cast<std::uint64_t>((i + 1) * base_step_size),
                 coverage_estimates[i] * bin_size);
}

// ADS: functions same, header different (above and this one)
static auto
write_predicted_coverage_curve(
  const std::string &outfile,                        //
  const double c_level,                              //
  const double base_step_size,                       //
  const std::uint32_t bin_size,                      //
  const std::vector<double> &cvrg_estimates,         //
  const std::vector<double> &cvrg_lower_ci_lognorm,  //
  const std::vector<double> &cvrg_upper_ci_lognorm) {
  static constexpr auto one_hundred = 100.0;

  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);

  const double percentile = one_hundred * c_level;
  // clang-format off
  std::println(out, "TOTAL_BASES\t"
               "EXPECTED_COVERED_BASES\t"
               "LOWER_{0}CI\t"
               "UPPER_{0}CI",
               percentile);
  // clang-format on

  std::println(out, "0\t0\t0\t0");
  for (auto i = 0U; i < std::size(cvrg_estimates); ++i)
    std::println(out, "{}\t{:.1f}\t{:.1f}\t{:.1f}",
                 static_cast<std::uint64_t>((i + 1) * base_step_size),
                 cvrg_estimates[i] * bin_size,
                 cvrg_lower_ci_lognorm[i] * bin_size,
                 cvrg_upper_ci_lognorm[i] * bin_size);
}

// NOLINTBEGIN(*-narrowing-conversions)

auto
gc_extrap_main(const std::span<char *> args) -> int {
  static constexpr auto cmd_name = "gc_extrap";
  static constexpr auto min_required_counts = 4;
  static constexpr auto max_iter_per_bootstrap = 10;
  static constexpr auto max_extrap_default = 1'000'000'000'000.0;

  const int argc = std::ssize(args);
  auto argv = std::data(args);

  std::string outfile;
  std::string infile;
  std::string histogram_outfile;

  // NOLINTBEGIN(*-avoid-magic-numbers)
  constexpr int diagonal{};
  std::uint32_t max_terms = 100;
  std::uint32_t bin_size = 10;
  double base_step_size = 1.0e8;
  std::uint32_t max_width = 10000;
  double max_extrap{max_extrap_default};
  std::uint32_t n_bootstraps = 100;
  std::uint32_t seed = 408;
  double c_level = 0.95;
  // NOLINTEND(*-avoid-magic-numbers)

  [[maybe_unused]] constexpr std::uint32_t n_threads{1};

  bool allow_defects{};
  bool verbose{};
  bool single_estimate{};

  CLI::App app{gc_extrap_about_msg};
  argv = app.ensure_utf8(argv);
  app.failure_message(
    [&](const CLI::App *a, const CLI::Error &b) -> std::string {
      return std::format("preseq command: {}\n", cmd_name) +
             CLI::FailureMessage::simple(a, b);
    });
  app.usage(std::format("Usage: preseq {} [OPTIONS]", cmd_name));
  if (argc >= 2)
    app.footer(gc_extrap_footer_msg);

  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  app.formatter(std::make_shared<preseq_formatter>());
  app.get_formatter()->column_width(24);
  app.get_formatter()->long_option_alignment_ratio(0.28);
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  // clang-format off
  app.set_help_flag("-h,--help", "print a detailed help message and exit");
  app.add_option("INPUT", infile, "input file (BED/BAM/SAM)")
    ->required()
    ->option_text(" ")
    ->check(CLI::ExistingFile)
    ->check(CLI::ReadPermissions);
  app.add_option("-o,--output", outfile, "coverage yield output file")
    ->required()
    ->option_text("FILE")
    ->check(CLI::WritePermissions);
  app.add_option("--hist-out", histogram_outfile, "write coverage counts histogram here")
    ->option_text("FILE");
  const auto max_width_opt =
    app.add_option("-w,--max-width", max_width, "max read length (longer are truncated)");
  app.add_option("--max_width", max_width, "for backwards compatibility")
    ->excludes(max_width_opt)
    ->group("Hidden");
  const auto bin_size_opt =
    app.add_option("-b,--bin-size", bin_size, "genomic bin size (partitions reads)");
  app.add_option("--bin_size", bin_size, "for backwards compatibility")
    ->excludes(bin_size_opt)
    ->group("Hidden");
  app.add_option("-e,--extrap", max_extrap, "maximum extrapolation (positive integer)")
    ->option_text("UINT")
    ->check(CLI::PositiveNumber);
  app.add_option("-s,--step", base_step_size, "step size in bases between extrapolations");
  app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
  app.add_option("-c,--cval", c_level, "level for confidence intervals");
  app.add_option("-x,--terms", max_terms, "maximum number of terms");
  app.add_option("-r,--seed", seed, "seed for random number generator");
  app.add_flag("-d,--defects", allow_defects,
               "make estimates without testing for defects");
  app.add_flag("-q,--quick", single_estimate, "no bootstraps for confidence intervals");
  app.add_flag("-v,--verbose", verbose, "print more info");
  // clang-format on

  if (argc < 2) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  if (!positive_integer(max_extrap)) {
    std::println("not integer: {}", max_extrap);
    return EXIT_FAILURE;
  }

  const auto bam_format_input = is_sam_or_bam_format(infile);
  const auto input_format = bam_format_input ? "BAM" : "BED";
  if (verbose)
    std::println("LOADING READS ({} format)", input_format);

  const auto [n_reads, coverage_hist] =
    bam_format_input
      ? load_coverage_counts_BAM(n_threads, infile, seed, bin_size, max_width)
      : load_coverage_counts(infile, seed, bin_size, max_width);

  const auto total_bins = get_counts_from_hist(coverage_hist);
  const auto distinct_bins =
    std::reduce(std::cbegin(coverage_hist), std::cend(coverage_hist));
  const double avg_bins_per_read = total_bins / n_reads;
  const double bin_step_size = base_step_size / bin_size;

  const std::uint32_t max_observed_count = std::size(coverage_hist) - 1;

  // ensure that the max terms are acceptable
  const auto first_zero = preseq::get_first_zero(coverage_hist);

  max_terms = std::min(max_terms, first_zero - 1);

  if (verbose)
    std::print("TOTAL READS         = {}\n"
               "BASE STEP SIZE      = {}\n"
               "BIN STEP SIZE       = {}\n"
               "TOTAL BINS          = {}\n"
               "BINS PER READ       = {}\n"
               "DISTINCT BINS       = {}\n"
               "TOTAL BASES         = {}\n"
               "TOTAL COVERED BASES = {}\n"
               "MAX COVERAGE COUNT  = {}\n"
               "COUNTS OF 1         = {}\n",  //
               n_reads,                       //
               base_step_size,                //
               bin_step_size,                 //
               total_bins,                    //
               avg_bins_per_read,             //
               distinct_bins,                 //
               total_bins * bin_size,         //
               distinct_bins * bin_size,      //
               max_observed_count,            //
               coverage_hist[1]               //
    );
  if (!histogram_outfile.empty())
    report_histogram(histogram_outfile, coverage_hist);
  // catch if all reads are distinct
  if (max_terms < min_required_counts)
    throw std::runtime_error("max count before zero is les than min required "
                             "count (4), sample not sufficiently deep or "
                             "duplicates removed");
  // check to make sure library is not overly saturated
  const double two_fold_extrap = preseq::good_toulmin_2x(coverage_hist);
  if (two_fold_extrap < 0.0)
    throw std::runtime_error("Library expected to saturate in doubling of "
                             "experiment size, unable to extrapolate");
  if (verbose)
    std::println("[ESTIMATING COVERAGE CURVE]");
  if (single_estimate) {
    const auto coverage_estimates = preseq::extrapolate_once(
      coverage_hist, max_terms, diagonal, allow_defects, bin_step_size,
      max_extrap / bin_size);
    if (coverage_estimates.empty()) {
      std::println("Single estimate failed. Run in full mode for estimates");
      return EXIT_FAILURE;
    }
    write_output(outfile, bin_size, base_step_size, coverage_estimates);
  }
  else {
    if (verbose)
      std::println("[BOOTSTRAPPING HISTOGRAM]");
    const std::uint32_t max_iter = max_iter_per_bootstrap * n_bootstraps;
    const auto cfa = preseq::cfa_t{
      .max_terms = max_terms,
      .diagonal = diagonal,
      .allow_defects = allow_defects,
      .rng_seed = seed,
      .target_bootstraps = n_bootstraps,
      .max_iterations = max_iter,
      .step_size = bin_step_size,
      .max_extrap = max_extrap / bin_size,
    };
    auto bootstrap_estimates =
      preseq::extrapolate_bootstrap(cfa, coverage_hist);
    if (verbose)
      std::println("[COMPUTING CONFIDENCE INTERVALS]");
    auto [medians, lower_ci_lognorms, upper_ci_lognorms] =
      median_and_ci_md(bootstrap_estimates, c_level);
    if (verbose)
      std::println("[WRITING OUTPUT]");
    write_predicted_coverage_curve(outfile, c_level, base_step_size, bin_size,
                                   medians, lower_ci_lognorms,
                                   upper_ci_lognorms);
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
