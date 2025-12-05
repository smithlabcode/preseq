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

#include "c_curve.hpp"

#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include "CLI11.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using std::accumulate;
using std::cbegin;
using std::cend;
using std::mt19937;
using std::size;
using std::size_t;
using std::string;
using std::uint32_t;
using std::vector;

template <typename T>
T
median_from_sorted_vector(const vector<T> &sorted_data, const size_t stride,
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
c_curve_main(int argc, char *argv[]) {
  try {
    bool verbose = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    uint32_t seed = 408;

    string outfile;
    string input_file_name;
    string histogram_outfile;

    double step_size = 1e6;
#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
    uint32_t n_threads{1};
#endif
    const auto description = R"(
Generate the complexity curve for data. This does not extrapolate, but instead
resamples from the given data.
)";
    CLI::App app{rlstrip(description)};
    argv = app.ensure_utf8(argv);
    app.usage("Usage: preseq c_curve [OPTIONS]");
    // if (argc >= 2)
    //   app.footer(description);

    // clang-format off
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "yield output file (default: stdout)");
    app.add_option("-s,--step", step_size, "step size in extrapolations");
    app.add_option("-P,--pe", PAIRED_END, "input is paired end read file");
    app.add_option("-H,--hist", HIST_INPUT, "input is a text file containing the observed histogram");
    app.add_option("-V,--vals", VALS_INPUT,
                   "input is a text file containing only the observed counts");
#ifdef HAVE_HTSLIB
    app.add_option("-B,--bam", BAM_FORMAT_INPUT, "input is in BAM format");
    app.add_option("-l,--seg_len", MAX_SEGMENT_LENGTH, "maximum segment length when merging paired end bam reads");
#endif
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_option("-v,--verbose", verbose, "print more info");
    // clang-format on

    if (argc < 3) {
      // std::println("{}", app.help());
      std::cout << app.help() << std::endl;
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    // Setup the random number generator
    mt19937 rng(seed);

    vector<double> counts_hist;
    size_t n_reads = 0;

    // LOAD VALUES
    if (HIST_INPUT) {
      if (verbose)
        std::cerr << "INPUT_HIST\n";
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (verbose)
        std::cerr << "VALS_INPUT\n";
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (verbose)
        std::cerr << "PAIRED_END_BAM_INPUT\n";
      n_reads = load_counts_BAM_pe(n_threads, input_file_name, counts_hist);
    }
    else if (BAM_FORMAT_INPUT) {
      if (verbose)
        std::cerr << "BAM_INPUT\n";
      n_reads = load_counts_BAM_se(n_threads, input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (verbose)
        std::cerr << "PAIRED_END_BED_INPUT\n";
      n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (verbose)
        std::cerr << "BED_INPUT\n";
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }

    const size_t max_observed_count = size(counts_hist) - 1;
    const double distinct_reads =
      accumulate(cbegin(counts_hist), cend(counts_hist), 0.0);

    const size_t total_reads = get_counts_from_hist(counts_hist);

    const size_t distinct_counts =
      std::count_if(cbegin(counts_hist), cend(counts_hist),
                    [](const double x) { return x > 0.0; });

    if (verbose)
      std::cerr << "TOTAL READS     = " << n_reads << '\n'
                << "COUNTS_SUM      = " << total_reads << '\n'
                << "DISTINCT READS  = " << distinct_reads << '\n'
                << "DISTINCT COUNTS = " << distinct_counts << '\n'
                << "MAX COUNT       = " << max_observed_count << '\n'
                << "COUNTS OF 1     = " << counts_hist[1] << '\n';

    if (!histogram_outfile.empty())
      report_histogram(histogram_outfile, counts_hist);

    const size_t upper_limit = n_reads;  // set upper limit equal to number of
                                         // molecules

    // setup for output of the complexity curve
    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    // prints the complexity curve
    out << "total_reads" << "\t" << "distinct_reads\n";
    out << 0 << '\t' << 0 << '\n';
    for (size_t i = step_size; i <= upper_limit; i += step_size) {
      if (verbose)
        std::cerr << "sample size: " << i << '\n';
      out << i << "\t"
          << interpolate_distinct(counts_hist, total_reads, distinct_reads, i)
          << '\n';
    }
  }
  catch (const std::exception &e) {
    std::cerr << "ERROR:\t" << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
