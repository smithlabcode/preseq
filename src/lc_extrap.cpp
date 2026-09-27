// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "lc_extrap.hpp"

#include "cli_common.hpp"
#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include <CLI11/CLI11.hpp>
#include <libpreseq.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iterator>
#include <memory>
#include <numeric>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

auto
lc_extrap_main(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
  try {
    static const std::size_t min_required_counts = 4;
    static const std::string min_required_counts_error_message =
      "max count before zero is less than min required count (" +
      std::to_string(min_required_counts) + ") duplicates removed";

    std::string outfile;
    std::string input_file_name;
    std::string histogram_outfile;

    // NOLINTBEGIN(*-avoid-magic-numbers)
    int diagonal{0};
    std::size_t orig_max_terms{100};
    std::uint64_t max_extrap_int{10'000'000'000};
    double step_size{1e6};
    std::size_t n_bootstraps{100};
    double c_level{0.95};
    std::uint32_t seed{408};
    // NOLINTEND(*-avoid-magic-numbers)

    // flags
    bool verbose{false};
    bool paired_end{false};
    bool single_estimate{false};
    bool allow_defects{false};

#ifdef HAVE_HTSLIB
    std::uint32_t n_threads{1};
#endif

    CLI::App app{lc_extrap_about_msg};
    argv = app.ensure_utf8(argv);
    app.formatter(std::make_shared<preseq_formatter>());
    app.usage("\nUsage: preseq lc_extrap [OPTIONS]");
    if (argc >= 3)
      app.footer(lc_extrap_footer_msg);

    // clang-format off
    app.set_help_flag("-h,--help", "print a detailed help message and exit");
    app.add_option("INPUT", input_file_name,
                   "input file name")
      ->required()
      ->option_text(" ")
      ->check(CLI::ExistingFile)
      // ->check(CLI::ReadPermission)
      ;
    app.add_option("-o,--output", outfile, "output filename (directory must exist)")
      ->option_text("FILE")
      // ->check(CLI::WritePermission)
      ->required();
    app.add_option("-e,--extrap", max_extrap_int, "maximum extrapolation (must be integer)")
      ->check(CLI::PositiveNumber);
    app.add_option("-s,--step", step_size, "extrapolation step size");
    app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--cval", c_level, "level for confidence intervals");
    app.add_option("-x,--terms", orig_max_terms, "maximum terms in estimator");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_flag("-p,--paired-end", paired_end, "input is paired end read file");
    app.add_flag("-q,--quick", single_estimate, "no bootstraps for confidence intervals");
    app.add_flag("-D,--defects", allow_defects, "no testing for defects");
    app.add_flag("-v,--verbose", verbose, "print more info");
    // clang-format on

    if (argc < 3) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const double max_extrap = max_extrap_int;

    const auto input_format = get_input_format_type(input_file_name);
    if (is_unknown(input_format)) {
      std::println("unknown input format");
      return EXIT_FAILURE;
    }

    if (verbose) {
      std::println("INPUT FORMAT: {}", to_string(input_format));
      if (is_bam(input_format) || is_bed(input_format))
        std::println("PAIRED END: {}", paired_end);
    }

    const auto [n_reads, counts_hist] = [&] {
      switch (input_format) {
      case input_format_type::hist:
        return load_histogram(input_file_name);
      case input_format_type::counts:
        return load_counts(input_file_name);
#ifdef HAVE_HTSLIB
      case input_format_type::bam:
        return paired_end ? load_counts_BAM_pe(n_threads, input_file_name)
                          : load_counts_BAM_se(n_threads, input_file_name);
#endif
      default:  // case input_format_type::vals:
        return paired_end ? load_counts_bed_pe(input_file_name)
                          : load_counts_bed_se(input_file_name);
      }
    }();

    const std::size_t max_observed_count = std::size(counts_hist) - 1;
    const double distinct_reads =
      std::reduce(std::cbegin(counts_hist), std::cend(counts_hist));

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::size_t first_zero = 1;
    while (first_zero < std::size(counts_hist) && counts_hist[first_zero] > 0)
      ++first_zero;

    // make sure the max terms is at most one less than the first zero
    orig_max_terms = std::min(orig_max_terms, first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

    const std::size_t distinct_counts =
      std::count_if(std::cbegin(counts_hist), std::cend(counts_hist),
                    [](const double x) { return x > 0.0; });

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

    if (!histogram_outfile.empty()) {
      report_histogram(histogram_outfile, counts_hist);
    }

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
      std::println("{}\t{}\t{}\t{}", orig_max_terms, diagonal, step_size,
                   max_extrap);
      const auto estimates =
        preseq::extrapolate_once(counts_hist, orig_max_terms, diagonal,
                                 allow_defects, step_size, max_extrap);
      if (estimates.empty())
        throw std::runtime_error("single estimate failed; run in full mode");

      std::ofstream out(outfile);
      if (!out)
        throw std::runtime_error("failed to open output file: " + outfile);

      std::println(out, "TOTAL_READS\tEXPECTED_DISTINCT");
      std::println(out, "0\t0");
      for (auto i = 0LU; i < std::size(estimates); ++i)
        std::println(out, "{}\t{}", (i + 1) * step_size, estimates[i]);
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
      auto bootstrap_estimates =
        preseq::extrapolate_bootstrap(cfa, counts_hist);
      if (verbose)
        std::println("[COMPUTING CONFIDENCE INTERVALS]");
      auto [medians, lower_ci_lognorms, upper_ci_lognorms] =
        median_and_ci_md(bootstrap_estimates, c_level);
      if (verbose)
        std::println("[WRITING OUTPUT]");
      write_predicted_complexity_curve(outfile, c_level, step_size, medians,
                                       lower_ci_lognorms, upper_ci_lognorms);
    }
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
