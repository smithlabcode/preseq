/* Copyright (C) 2013-2025 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Authors: Timothy Daley and Andrew Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

static constexpr auto about_msg = R"(
Generate the full observed complexity curve for data. This does not
extrapolate, but instead resamples from the given data.
)";

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
#include <string>
#include <vector>

int
c_curve_main(int argc, char *argv[]) {
  try {
    bool verbose = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    std::uint32_t seed = 408;

    std::string outfile;
    std::string input_file_name;
    std::string histogram_outfile;

    double step_size = 1e6;
#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    std::size_t MAX_SEGMENT_LENGTH = 5000;
    std::uint32_t n_threads{1};
#endif
    CLI::App app{rlstrip(about_msg)};
    argv = app.ensure_utf8(argv);
    app.usage("\nUsage: preseq c_curve [OPTIONS]");
    // if (argc >= 3)
    //   app.footer(description);

    // clang-format off
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "yield output file (default: stdout)");
    app.add_option("-s,--step", step_size, "step size in extrapolations");
    app.add_flag("-P,--pe", PAIRED_END, "input is paired end read file");
    app.add_flag("-H,--hist", HIST_INPUT,
                 "input is a text file containing the observed histogram");
    app.add_flag("-V,--vals", VALS_INPUT,
                   "input is a text file containing only the observed counts");
#ifdef HAVE_HTSLIB
    app.add_flag("-B,--bam", BAM_FORMAT_INPUT, "input is in BAM format");
    app.add_option("-l,--seg_len", MAX_SEGMENT_LENGTH,
                   "maximum segment length when merging paired end bam reads");
#endif
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_flag("-v,--verbose", verbose, "print more info");
    // clang-format on

    if (argc < 3) {
      // std::println("{}", app.help());
      std::cout << app.help() << std::endl;
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    std::vector<double> counts_hist;
    std::size_t n_reads = 0;

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

    const auto max_observed_count = std::size(counts_hist) - 1;
    const auto distinct_reads =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    const auto total_reads = get_counts_from_hist(counts_hist);
    const auto distinct_counts =
      std::count_if(std::cbegin(counts_hist), std::cend(counts_hist),
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

    // set upper limit equal to number of molecules
    const std::size_t upper_limit = n_reads;

    // setup for output of the complexity curve
    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    // prints the complexity curve
    out << "total_reads" << "\t" << "distinct_reads\n";
    out << 0 << '\t' << 0 << '\n';
    for (std::size_t i = step_size; i <= upper_limit; i += step_size) {
      if (verbose)
        std::cerr << "sample size: " << i << '\n';
      out << i << '\t'
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
