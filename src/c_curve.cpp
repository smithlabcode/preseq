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

#include "c_curve.hpp"

#include "common.hpp"
#include "continued_fraction.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

#include <GenomicRegion.hpp>
#include <OptionParser.hpp>
#include <smithlab_os.hpp>
#include <smithlab_utils.hpp>

#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using std::array;
using std::cerr;
using std::endl;
using std::isfinite;
using std::max;
using std::min;
using std::mt19937;
using std::runtime_error;
using std::setprecision;
using std::size;
using std::string;
using std::to_string;
using std::uint64_t;
using std::unordered_map;
using std::vector;

template <typename T>
T
median_from_sorted_vector(const vector<T> sorted_data, const size_t stride,
                          const size_t n) {
  if (n == 0 || sorted_data.empty())
    return 0.0;

  const size_t lhs = (n - 1) / 2;
  const size_t rhs = n / 2;

  if (lhs == rhs)
    return sorted_data[lhs * stride];

  return (sorted_data[lhs * stride] + sorted_data[rhs * stride]) / 2.0;
}

int
c_curve_main(const int argc, const char *argv[]) {
  try {
    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    uint64_t seed = 408;

    string outfile;

    size_t upper_limit = 0;
    double step_size = 1e6;
#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
    uint32_t n_threads{1};
#endif

    const string description = R"(
Generate the complexity curve for data. This does not extrapolate, \
but instead resamples from the given data.)";

    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]), description, "<input-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("step", 's', "step size in extrapolations", false,
                      step_size);
    opt_parse.add_opt("verbose", 'v', "print more information", false, VERBOSE);
    opt_parse.add_opt("pe", 'P', "input is paired end read file", false,
                      PAIRED_END);
    opt_parse.add_opt("hist", 'H',
                      "input is a text file containing the observed histogram",
                      false, HIST_INPUT);
    opt_parse.add_opt(
      "vals", 'V', "input is a text file containing only the observed counts",
      false, VALS_INPUT);
#ifdef HAVE_HTSLIB
    opt_parse.add_opt("bam", 'B', "input is in BAM format", false,
                      BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l',
                      "maximum segment length when merging "
                      "paired end bam reads",
                      false, MAX_SEGMENT_LENGTH);
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
    const string input_file_name = leftover_args.front();
    /******************************************************************/

    // Setup the random number generator
    srand(time(0) + getpid());  // give the random fxn a new seed
    mt19937 rng(seed);

    vector<double> counts_hist;
    size_t n_reads = 0;

    // LOAD VALUES
    if (HIST_INPUT) {
      if (VERBOSE)
        cerr << "INPUT_HIST" << endl;
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
      n_reads = load_counts_BAM_pe(n_threads, input_file_name, counts_hist);
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        cerr << "BAM_INPUT" << endl;
      n_reads = load_counts_BAM_se(n_threads, input_file_name, counts_hist);
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

    const size_t max_observed_count = counts_hist.size() - 1;
    const double distinct_reads =
      accumulate(begin(counts_hist), end(counts_hist), 0.0);

    const size_t total_reads = get_counts_from_hist(counts_hist);

    const size_t distinct_counts =
      std::count_if(begin(counts_hist), end(counts_hist),
                    [](const double x) { return x > 0.0; });

    if (VERBOSE)
      cerr << "TOTAL READS     = " << n_reads << endl
           << "COUNTS_SUM      = " << total_reads << endl
           << "DISTINCT READS  = " << distinct_reads << endl
           << "DISTINCT COUNTS = " << distinct_counts << endl
           << "MAX COUNT       = " << max_observed_count << endl
           << "COUNTS OF 1     = " << counts_hist[1] << endl;

    if (VERBOSE) {
      // output the original histogram
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0)
          cerr << i << '\t' << static_cast<size_t>(counts_hist[i]) << endl;
      cerr << endl;
    }

    if (upper_limit == 0)
      upper_limit = n_reads;  // set upper limit to equal the number of
                              // molecules

    // handles output of c_curve
    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    // prints the complexity curve
    out << "total_reads" << "\t" << "distinct_reads" << endl;
    out << 0 << '\t' << 0 << endl;
    for (size_t i = step_size; i <= upper_limit; i += step_size) {
      if (VERBOSE)
        cerr << "sample size: " << i << endl;
      out << i << "\t"
          << interpolate_distinct(counts_hist, total_reads, distinct_reads, i)
          << endl;
    }
  }
  catch (runtime_error &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
