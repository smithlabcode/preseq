/* Copyright (C) 2013-2026 Andrew D. Smith and Timothy Daley
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
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
#include <iterator>
#include <numeric>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

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

auto
gc_extrap_main(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
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

    CLI::App app{rlstrip(gc_extrap_about_msg)};
    argv = app.ensure_utf8(argv);
    app.usage("\nUsage: preseq gc_extrap [OPTIONS]");
    if (argc >= 3)
      app.footer(rlstrip(gc_extrap_footer_msg));

    // clang-format off
    app.set_help_flag("-h,--help", "print a detailed help message and exit");
    app.add_option("INPUT", infile, "input file name")
      ->required()
      ->option_text(" ")
      ->check(CLI::ExistingFile)
      // ->check(CLI::ReadPermission)
      ;
    app.add_option("-o,--output", outfile, "coverage yield output file")
      ->option_text("FILE")
      // ->check(CLI::WritePermission)
      ->required();
    app.add_option("--hist-out", histogram_outfile, "output histogram to this file")
      ->option_text("FILE");
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
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const double max_extrap = max_extrap_int;

    const auto bam_format_input = is_sam_or_bam_format(infile);

    const auto input_format = bam_format_input ? "BAM" : "BED";
    if (verbose)
      std::println("LOADING READS ({} format)", input_format);

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
      std::println("TOTAL READS         = {}"
                   "BASE STEP SIZE      = {}"
                   "BIN STEP SIZE       = {}"
                   "TOTAL BINS          = {}"
                   "BINS PER READ       = {}"
                   "DISTINCT BINS       = {}"
                   "TOTAL BASES         = {}"
                   "TOTAL COVERED BASES = {}"
                   "MAX COVERAGE COUNT  = {}"
                   "COUNTS OF 1         = {}",  //
                   n_reads,                     //
                   base_step_size,              //
                   bin_step_size,               //
                   total_bins,                  //
                   avg_bins_per_read,           //
                   distinct_bins,               //
                   total_bins * bin_size,       //
                   distinct_bins * bin_size,    //
                   max_observed_count,          //
                   coverage_hist[1]             //
      );

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
      std::println("[ESTIMATING COVERAGE CURVE]");

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
        std::println("[BOOTSTRAPPING HISTOGRAM]");

      const std::uint32_t max_iter = max_iter_per_bootstrap * n_bootstraps;

      std::vector<std::vector<double>> bootstrap_estimates;
      extrap_bootstrap(verbose, allow_defects, seed, coverage_hist,
                       n_bootstraps, orig_max_terms, diagonal, bin_step_size,
                       max_extrap / bin_size, max_iter, bootstrap_estimates);

      if (verbose)
        std::println("[COMPUTING CONFIDENCE INTERVALS]");
      std::vector<double> coverage_upper_ci_lognorm;
      std::vector<double> coverage_lower_ci_lognorm;
      vector_median_and_ci(bootstrap_estimates, c_level, coverage_estimates,
                           coverage_lower_ci_lognorm,
                           coverage_upper_ci_lognorm);

      if (verbose)
        std::println("[WRITING OUTPUT]");

      write_predicted_coverage_curve(
        outfile, c_level, base_step_size, bin_size, coverage_estimates,
        coverage_lower_ci_lognorm, coverage_upper_ci_lognorm);
    }
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
