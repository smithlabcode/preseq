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

#include "lc_extrap.hpp"

#include "common.hpp"

#include "load_data_for_complexity.hpp"

#include "OptionParser.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using std::begin;
using std::cerr;
using std::end;
using std::endl;
using std::runtime_error;
using std::size_t;
using std::string;
using std::to_string;
using std::vector;

namespace fs = std::filesystem;

int
lc_extrap_main(const int argc, const char **argv) {
  try {
    static const size_t min_required_counts = 4;
    static const string min_required_counts_error_message =
      "max count before zero is less than min required count (" +
      to_string(min_required_counts) + ") duplicates removed";

    string outfile;

    size_t orig_max_terms = 100;
    double max_extrap = 1.0e10;
    double step_size = 1e6;
    size_t n_bootstraps = 100;
    int diagonal = 0;
    double c_level = 0.95;
    uint64_t seed = 408;

    /* FLAGS */
    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool SINGLE_ESTIMATE = false;
    bool allow_defects = false;

#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    string description =
      R"(Extrapolate the complexity of a library. This is the approach   \
described in Daley & Smith (2013). The method applies rational  \
function approximation via continued fractions with the         \
original goal of estimating the number of distinct reads that a \
sequencing library would yield upon deeper sequencing. This     \
method has been used for many different purposes since then.
)";

    /********** GET COMMAND LINE ARGUMENTS  FOR LC EXTRAP ***********/

    OptionParser opt_parse(fs::path(argv[1]).filename(), description,
                           "<input-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("extrap", 'e', "maximum extrapolation", false,
                      max_extrap);
    opt_parse.add_opt("step", 's', "extrapolation step size", false, step_size);
    opt_parse.add_opt("boots", 'n', "number of bootstraps", false,
                      n_bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals", false,
                      c_level);
    opt_parse.add_opt("terms", 'x', "maximum terms in estimator", false,
                      orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more info", false, VERBOSE);
#ifdef HAVE_HTSLIB
    opt_parse.add_opt("bam", 'B', "input is in BAM format", false,
                      BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l',
                      "maximum segment length when merging "
                      "paired end bam reads",
                      false, MAX_SEGMENT_LENGTH);
#endif
    opt_parse.add_opt("pe", 'P', "input is paired end read file", false,
                      PAIRED_END);
    opt_parse.add_opt(
      "vals", 'V', "input is a text file containing only the observed counts",
      false, VALS_INPUT);
    opt_parse.add_opt("hist", 'H',
                      "input is a text file containing the observed histogram",
                      false, HIST_INPUT);
    opt_parse.add_opt("quick", 'Q',
                      "quick mode (no bootstraps) for confidence intervals",
                      false, SINGLE_ESTIMATE);
    opt_parse.add_opt("defects", 'D', "no testing for defects", false,
                      allow_defects);
    opt_parse.add_opt("seed", 'r', "seed for random number generator", false,
                      seed);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    // ADS: suspect bug below; "-about" isn't working.
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
    const string input_file_name = leftover_args.front();
    /******************************************************************/

    vector<double> counts_hist;
    size_t n_reads = 0;

    /************ loading input ***************************************/
    if (HIST_INPUT) {
      if (VERBOSE)
        cerr << "HIST_INPUT" << endl;
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (VERBOSE)
        cerr << "VALS_INPUT" << endl;
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (VERBOSE)
        cerr << "PAIRED_END_BAM_INPUT" << endl;
      const size_t MAX_READS_TO_HOLD = 5000000;
      size_t n_paired = 0;
      size_t n_mates = 0;
      n_reads =
        load_counts_BAM_pe(VERBOSE, input_file_name, MAX_SEGMENT_LENGTH,
                           MAX_READS_TO_HOLD, n_paired, n_mates, counts_hist);
      if (VERBOSE) {
        cerr << "MERGED PAIRED END READS = " << n_paired << endl;
        cerr << "MATES PROCESSED = " << n_mates << endl;
      }
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        cerr << "BAM_INPUT" << endl;
      n_reads = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (VERBOSE)
        cerr << "PAIRED_END_BED_INPUT" << endl;
      n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }
    /************ done loading input **********************************/

    const size_t max_observed_count = counts_hist.size() - 1;
    const double distinct_reads =
      std::accumulate(begin(counts_hist), end(counts_hist), 0.0);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t first_zero = 1;
    while (first_zero < counts_hist.size() && counts_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = std::min(orig_max_terms, first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

    const size_t distinct_counts =
      std::count_if(begin(counts_hist), end(counts_hist),
                    [](const double x) { return x > 0.0; });

    if (VERBOSE)
      cerr << "TOTAL READS     = " << n_reads << endl
           << "DISTINCT READS  = " << distinct_reads << endl
           << "DISTINCT COUNTS = " << distinct_counts << endl
           << "MAX COUNT       = " << max_observed_count << endl
           << "COUNTS OF 1     = " << counts_hist[1] << endl
           << "MAX TERMS       = " << orig_max_terms << endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0)
          cerr << i << '\t' << static_cast<size_t>(counts_hist[i]) << endl;
      cerr << endl;
    }

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(counts_hist);
    if (two_fold_extrap < 0.0)
      throw runtime_error("Saturation expected at double initial sample size. "
                          "Unable to extrapolate.");

    // const size_t total_reads = get_counts_from_hist(counts_hist);

    // assert(total_reads == n_reads); // ADS: why commented out?

    // check that min required count is satisfied
    if (orig_max_terms < min_required_counts)
      throw runtime_error(min_required_counts_error_message);

    if (VERBOSE)
      cerr << "[ESTIMATING YIELD CURVE]" << endl;
    vector<double> yield_estimates;

    if (SINGLE_ESTIMATE) {
      const bool single_estimate_success = extrap_single_estimate(
        VERBOSE, allow_defects, counts_hist, orig_max_terms, diagonal,
        step_size, max_extrap, yield_estimates);
      // IF FAILURE, EXIT
      if (!single_estimate_success)
        throw runtime_error("single estimate failed, run "
                            "full mode for estimates");

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS\tEXPECTED_DISTINCT" << endl;
      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << endl;
      for (size_t i = 0; i < yield_estimates.size(); ++i)
        out << (i + 1) * step_size << '\t' << yield_estimates[i] << endl;
    }
    else {
      if (VERBOSE)
        cerr << "[BOOTSTRAPPING HISTOGRAM]" << endl;

      const size_t max_iter = 100 * n_bootstraps;

      vector<vector<double>> bootstrap_estimates;
      extrap_bootstrap(VERBOSE, allow_defects, seed, counts_hist, n_bootstraps,
                       orig_max_terms, diagonal, step_size, max_extrap,
                       max_iter, bootstrap_estimates);

      if (VERBOSE)
        cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      // yield ci
      vector<double> yield_upper_ci_lognorm, yield_lower_ci_lognorm;

      vector_median_and_ci(bootstrap_estimates, c_level, yield_estimates,
                           yield_lower_ci_lognorm, yield_upper_ci_lognorm);

      /////////////////////////////////////////////////////////////////////
      if (VERBOSE)
        cerr << "[WRITING OUTPUT]" << endl;

      write_predicted_complexity_curve(outfile, c_level, step_size,
                                       yield_estimates, yield_lower_ci_lognorm,
                                       yield_upper_ci_lognorm);
    }
  }
  catch (std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
