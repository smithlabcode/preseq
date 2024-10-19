/* Copyright (C) 2013-2024 University of Southern California and
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

#include <OptionParser.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using std::cbegin;
using std::cend;
using std::cerr;
using std::endl;
using std::min;
using std::runtime_error;
using std::string;
using std::vector;

namespace fs = std::filesystem;

// ADS: functions same, header different (above and this one)
static void
write_predicted_coverage_curve(const string &outfile, const double c_level,
                               const double base_step_size,
                               const size_t bin_size,
                               const vector<double> &cvrg_estimates,
                               const vector<double> &cvrg_lower_ci_lognorm,
                               const vector<double> &cvrg_upper_ci_lognorm) {
  static constexpr double one_hundred = 100.0;
  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  const double percentile = one_hundred * c_level;
  // clang-format off
  out << "TOTAL_BASES" << '\t'
      << "EXPECTED_COVERED_BASES" << '\t'
      << "LOWER_" << percentile << "%CI" << '\t'
      << "UPPER_" << percentile << "%CI"
      << endl;
  // clang-format on

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < cvrg_estimates.size(); ++i)
    out << (i + 1) * base_step_size << '\t' << cvrg_estimates[i] * bin_size
        << '\t' << cvrg_lower_ci_lognorm[i] * bin_size << '\t'
        << cvrg_upper_ci_lognorm[i] * bin_size << endl;
}

int
gc_extrap_main(const int argc, const char *argv[]) {
  try {
    const size_t MIN_REQUIRED_COUNTS = 4;

    string outfile;
    string histogram_outfile;

    int diagonal = 0;
    size_t orig_max_terms = 100;
    size_t bin_size = 10;
    bool verbose = false;
    double base_step_size = 1.0e8;
    size_t max_width = 10000;
    bool SINGLE_ESTIMATE = false;
    double max_extrap = 1.0e12;
    size_t n_bootstraps = 100;
    uint64_t seed = 408;
    bool allow_defects = false;

    bool NO_SEQUENCE = false;
    double c_level = 0.95;
#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    uint32_t n_threads{1};
#endif

    const string description = R"(
Extrapolate the size of the covered genome by mapped reads. This
approach is described in Daley & Smith (2014). The method is the same
as for lc_extrap: using rational function approximation to a
power-series expansion for the number of "unobserved" bases in the
initial sample. The gc_extrap method is adapted to deal with
individual nucleotides rather than distinct reads.
)";
    string program_name = fs::path(argv[0]).filename();
    program_name += " " + string(argv[1]);

    // ********* GET COMMAND LINE ARGUMENTS  FOR GC EXTRAP **********
    OptionParser opt_parse(program_name, description, "<input-file>");
    opt_parse.add_opt("output", 'o',
                      "coverage yield output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("max_width", 'w',
                      "max fragment length, "
                      "set equal to read length for single end reads",
                      false, max_width);
    opt_parse.add_opt("bin_size", 'b', "bin size", false, bin_size);
    opt_parse.add_opt("extrap", 'e', "maximum extrapolation in base pairs",
                      false, max_extrap);
    opt_parse.add_opt("step", 's', "step size in bases between extrapolations",
                      false, base_step_size);
    opt_parse.add_opt("bootstraps", 'n', "number of bootstraps", false,
                      n_bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals", false,
                      c_level);
    opt_parse.add_opt("terms", 'x', "maximum number of terms", false,
                      orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", false, verbose);
    opt_parse.add_opt("hist-out", '\0', "output histogram to this file", false,
                      histogram_outfile);
    opt_parse.add_opt("bed", 'B',
                      "input is in bed format without sequence information",
                      false, NO_SEQUENCE);
    opt_parse.add_opt("quick", 'Q',
                      "quick mode: run gc_extrap without "
                      "bootstrapping for confidence intervals",
                      false, SINGLE_ESTIMATE);
    opt_parse.add_opt("defects", 'D',
                      "defects mode to extrapolate without testing for defects",
                      false, allow_defects);
#ifdef HAVE_HTSLIB
    opt_parse.add_opt("bam", '\0', "input is in BAM format", false,
                      BAM_FORMAT_INPUT);
    opt_parse.add_opt("threads", 't', "number of threads for decompressing BAM",
                      false, n_threads);
#endif
    opt_parse.add_opt("seed", 'r', "seed for random number generator", false,
                      seed);
    opt_parse.set_show_defaults();

    vector<string> leftover_args;
    opt_parse.parse(argc - 1, argv + 1, leftover_args);
    if (argc == 2 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string infile = leftover_args.front();
    // ****************************************************************

    vector<double> coverage_hist;
    size_t n_reads = 0;
    if (verbose)
      cerr << "LOADING READS" << endl;

    if (NO_SEQUENCE) {
      if (verbose)
        cerr << "BED FORMAT" << endl;
      n_reads = load_coverage_counts_GR(infile, seed, bin_size, max_width,
                                        coverage_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT) {
      if (verbose)
        cerr << "BAM_INPUT" << endl;
      n_reads = load_coverage_counts_BAM(n_threads, infile, seed, bin_size,
                                         max_width, coverage_hist);
    }
#endif
    else {
      if (verbose)
        cerr << "MAPPED READ FORMAT" << endl;
      n_reads = load_coverage_counts_MR(infile, seed, bin_size, max_width,
                                        coverage_hist);
    }

    const double total_bins = get_counts_from_hist(coverage_hist);

    if (verbose)
      report_histogram(histogram_outfile, coverage_hist);

    const double distinct_bins =
      accumulate(cbegin(coverage_hist), cend(coverage_hist), 0.0);

    const double avg_bins_per_read = total_bins / n_reads;
    const double bin_step_size = base_step_size / bin_size;

    const size_t max_observed_count = coverage_hist.size() - 1;

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t first_zero = 1;
    while (first_zero < coverage_hist.size() && coverage_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = min(orig_max_terms, first_zero - 1);

    if (verbose)
      cerr << "TOTAL READS         = " << n_reads << endl
           << "BASE STEP SIZE      = " << base_step_size << endl
           << "BIN STEP SIZE       = " << bin_step_size << endl
           << "TOTAL BINS          = " << total_bins << endl
           << "BINS PER READ       = " << avg_bins_per_read << endl
           << "DISTINCT BINS       = " << distinct_bins << endl
           << "TOTAL BASES         = " << total_bins * bin_size << endl
           << "TOTAL COVERED BASES = " << distinct_bins * bin_size << endl
           << "MAX COVERAGE COUNT  = " << max_observed_count << endl
           << "COUNTS OF 1         = " << coverage_hist[1] << endl;

    if (verbose) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED BIN COUNTS (" << coverage_hist.size() << ")" << endl;
      for (size_t i = 0; i < coverage_hist.size(); i++)
        if (coverage_hist[i] > 0)
          cerr << i << '\t' << coverage_hist[i] << endl;
      cerr << endl;
    }

    // catch if all reads are distinct
    if (orig_max_terms < MIN_REQUIRED_COUNTS)
      throw runtime_error("max count before zero is les than min required "
                          "count (4), sample not sufficiently deep or "
                          "duplicates removed");

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(coverage_hist);
    if (two_fold_extrap < 0.0)
      throw runtime_error("Library expected to saturate in doubling of "
                          "experiment size, unable to extrapolate");

    if (verbose)
      cerr << "[ESTIMATING COVERAGE CURVE]" << endl;

    vector<double> coverage_estimates;

    if (SINGLE_ESTIMATE) {
      bool SINGLE_ESTIMATE_SUCCESS = extrap_single_estimate(
        verbose, allow_defects, coverage_hist, orig_max_terms, diagonal,
        bin_step_size, max_extrap / bin_size, coverage_estimates);
      // IF FAILURE, EXIT
      if (!SINGLE_ESTIMATE_SUCCESS)
        throw runtime_error("SINGLE ESTIMATE FAILED, NEED TO RUN IN "
                            "FULL MODE FOR ESTIMATES");

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile);
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_BASES\tEXPECTED_DISTINCT" << endl;

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << endl;
      for (size_t i = 0; i < coverage_estimates.size(); ++i)
        out << (i + 1) * base_step_size << '\t'
            << coverage_estimates[i] * bin_size << endl;
    }
    else {
      if (verbose)
        cerr << "[BOOTSTRAPPING HISTOGRAM]" << endl;

      const size_t max_iter = 10 * n_bootstraps;

      vector<vector<double>> bootstrap_estimates;
      extrap_bootstrap(verbose, allow_defects, seed, coverage_hist,
                       n_bootstraps, orig_max_terms, diagonal, bin_step_size,
                       max_extrap / bin_size, max_iter, bootstrap_estimates);

      if (verbose)
        cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      vector<double> coverage_upper_ci_lognorm, coverage_lower_ci_lognorm;
      vector_median_and_ci(bootstrap_estimates, c_level, coverage_estimates,
                           coverage_lower_ci_lognorm,
                           coverage_upper_ci_lognorm);

      if (verbose)
        cerr << "[WRITING OUTPUT]" << endl;

      write_predicted_coverage_curve(
        outfile, c_level, base_step_size, bin_size, coverage_estimates,
        coverage_lower_ci_lognorm, coverage_upper_ci_lognorm);
    }
  }
  catch (std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
