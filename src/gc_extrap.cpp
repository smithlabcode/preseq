/* Copyright (C) 2013-2025 University of Southern California and
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

#include "gc_extrap.hpp"

#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include "CLI11/CLI11.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

static auto
write_output(const std::string &outfile, const std::uint32_t bin_size,
             const std::uint32_t base_step_size,
             const std::vector<double> &coverage_estimates) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);

  out << "TOTAL_BASES\tEXPECTED_DISTINCT\n";

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out << 0 << '\t' << 0 << '\n';
  for (auto i = 0u; i < std::size(coverage_estimates); ++i)
    out << (i + 1) * base_step_size << '\t' << coverage_estimates[i] * bin_size
        << '\n';
}

// ADS: functions same, header different (above and this one)
static auto
write_predicted_coverage_curve(
  const std::string &outfile, const double c_level, const double base_step_size,
  const std::uint32_t bin_size, const std::vector<double> &cvrg_estimates,
  const std::vector<double> &cvrg_lower_ci_lognorm,
  const std::vector<double> &cvrg_upper_ci_lognorm) {
  static constexpr double one_hundred = 100.0;
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);

  const double percentile = one_hundred * c_level;
  // clang-format off
  out << "TOTAL_BASES" << '\t'
      << "EXPECTED_COVERED_BASES" << '\t'
      << "LOWER_" << percentile << "_CI" << '\t'
      << "UPPER_" << percentile << "_CI\n";
  // clang-format on

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\n';
  for (auto i = 0u; i < std::size(cvrg_estimates); ++i) {
    // clang-format off
    out << (i + 1) * base_step_size << '\t'
        << cvrg_estimates[i] * bin_size << '\t'
        << cvrg_lower_ci_lognorm[i] * bin_size << '\t'
        << cvrg_upper_ci_lognorm[i] * bin_size << '\n';
    // clang-format on
  }
}

auto
gc_extrap::main(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
  try {
    static constexpr auto min_required_counts = 4;
    static constexpr auto max_iter_per_bootstrap = 10;

    std::string outfile;
    std::string infile;
    std::string histogram_outfile;

    // NOLINTBEGIN(*-avoid-magic-numbers)
    int diagonal = 0;
    std::uint32_t orig_max_terms = 100;
    std::uint32_t bin_size = 10;
    double base_step_size = 1.0e8;
    std::uint32_t max_width = 10000;
    std::uint64_t max_extrap_int{1'000'000'000'000};
    std::uint32_t n_bootstraps = 100;
    std::uint32_t seed = 408;
    double c_level = 0.95;
    // NOLINTEND(*-avoid-magic-numbers)

#ifdef HAVE_HTSLIB
    std::uint32_t n_threads{1};
#endif

    bool allow_defects{false};
    bool verbose{false};
    bool single_estimate{false};

    CLI::App app{rlstrip(about_msg)};
    argv = app.ensure_utf8(argv);
    app.usage("\nUsage: preseq gc_extrap [OPTIONS]");
    if (argc >= 3)
      app.footer(rlstrip(footer_msg));

    // clang-format off
    app.set_help_flag("-h,--help", "print a detailed help message and exit");
    app.add_option("-i,--input", infile, "input file")
      ->option_text("FILE")
      ->required()
      // ->check(CLI::ReadPermission)
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "coverage yield output file")
      ->option_text("FILE")
      // ->check(CLI::WritePermission)
      ->required();
    app.add_option("-w,--max_width", max_width,
                   "max fragment length, set equal to read length for single end reads");
    app.add_option("-b,--bin_size", bin_size, "bin size");
    app.add_option("-e,--extrap", max_extrap_int, "maximum extrapolation (must be integer)")
      ->check(CLI::PositiveNumber);
    app.add_option("-s,--step", base_step_size, "step size in bases between extrapolations");
    app.add_option("-n,--bootstraps", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--cval", c_level, "level for confidence intervals");
    app.add_option("-x,--terms", orig_max_terms, "maximum number of terms");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_flag("-q,--quick", single_estimate, "no bootstraps for confidence intervals");
    app.add_flag("-D,--defects", allow_defects,
                 "defects mode to extrapolate without testing for defects");
    app.add_flag("-v,--verbose", verbose, "print more info");
    // clang-format on

    if (argc < 3) {
      // std::println("{}", app.help());
      std::cout << app.help() << '\n';
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const double max_extrap = max_extrap_int;

    const auto bam_format_input = is_sam_or_bam_format(infile);

    const auto input_format = bam_format_input ? "BAM" : "BED";
    if (verbose)
      std::cerr << "LOADING READS (" << input_format << " format)\n";

    const auto [n_reads, coverage_hist] =
#ifdef HAVE_HTSLIB
      bam_format_input
        ? load_coverage_counts_BAM(n_threads, infile, seed, bin_size, max_width)
        :
#endif
        load_coverage_counts(infile, seed, bin_size, max_width);

    const auto total_bins = get_counts_from_hist(coverage_hist);
    const auto distinct_bins = std::accumulate(std::cbegin(coverage_hist),
                                               std::cend(coverage_hist), 0.0);
    const double avg_bins_per_read = total_bins / n_reads;
    const double bin_step_size = base_step_size / bin_size;

    const std::uint32_t max_observed_count = std::size(coverage_hist) - 1;

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::uint32_t first_zero{1};
    while (first_zero < std::size(coverage_hist) &&
           coverage_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = std::min(orig_max_terms, first_zero - 1);

    if (verbose)
      std::cerr << "TOTAL READS         = " << n_reads << '\n'
                << "BASE STEP SIZE      = " << base_step_size << '\n'
                << "BIN STEP SIZE       = " << bin_step_size << '\n'
                << "TOTAL BINS          = " << total_bins << '\n'
                << "BINS PER READ       = " << avg_bins_per_read << '\n'
                << "DISTINCT BINS       = " << distinct_bins << '\n'
                << "TOTAL BASES         = " << total_bins * bin_size << '\n'
                << "TOTAL COVERED BASES = " << distinct_bins * bin_size << '\n'
                << "MAX COVERAGE COUNT  = " << max_observed_count << '\n'
                << "COUNTS OF 1         = " << coverage_hist[1] << '\n';

    if (!histogram_outfile.empty())
      report_histogram(histogram_outfile, coverage_hist);

    // catch if all reads are distinct
    if (orig_max_terms < min_required_counts)
      throw std::runtime_error("max count before zero is les than min required "
                               "count (4), sample not sufficiently deep or "
                               "duplicates removed");

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(coverage_hist);
    if (two_fold_extrap < 0.0)
      throw std::runtime_error("Library expected to saturate in doubling of "
                               "experiment size, unable to extrapolate");

    if (verbose)
      std::cerr << "[ESTIMATING COVERAGE CURVE]\n";

    std::vector<double> coverage_estimates;

    if (single_estimate) {
      const auto success = extrap_single_estimate(
        verbose, allow_defects, coverage_hist, orig_max_terms, diagonal,
        bin_step_size, max_extrap / bin_size, coverage_estimates);
      // IF FAILURE, EXIT
      if (!success)
        throw std::runtime_error(
          "Single estimate failed. Run in full mode for estimates");

      write_output(outfile, bin_size, base_step_size, coverage_estimates);
    }
    else {
      if (verbose)
        std::cerr << "[BOOTSTRAPPING HISTOGRAM]\n";

      const std::uint32_t max_iter = max_iter_per_bootstrap * n_bootstraps;

      std::vector<std::vector<double>> bootstrap_estimates;
      extrap_bootstrap(verbose, allow_defects, seed, coverage_hist,
                       n_bootstraps, orig_max_terms, diagonal, bin_step_size,
                       max_extrap / bin_size, max_iter, bootstrap_estimates);

      if (verbose)
        std::cerr << "[COMPUTING CONFIDENCE INTERVALS]\n";
      std::vector<double> coverage_upper_ci_lognorm, coverage_lower_ci_lognorm;
      vector_median_and_ci(bootstrap_estimates, c_level, coverage_estimates,
                           coverage_lower_ci_lognorm,
                           coverage_upper_ci_lognorm);

      if (verbose)
        std::cerr << "[WRITING OUTPUT]\n";

      write_predicted_coverage_curve(
        outfile, c_level, base_step_size, bin_size, coverage_estimates,
        coverage_lower_ci_lognorm, coverage_upper_ci_lognorm);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
