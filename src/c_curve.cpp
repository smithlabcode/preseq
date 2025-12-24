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

#include "c_curve.hpp"

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

auto
c_curve::main(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
  try {
    std::uint32_t seed = 408;  // NOLINT(*-avoid-magic-numbers)
    double step_size = 1e6;    // NOLINT(*-avoid-magic-numbers)

    std::string outfile;
    std::string infile;
    std::string histogram_outfile;

    bool verbose = false;
    bool paired_end = false;

#ifdef HAVE_HTSLIB
    std::uint32_t n_threads{1};
#endif

    CLI::App app{rlstrip(about_msg)};
    argv = app.ensure_utf8(argv);
    app.usage("\nUsage: preseq c_curve [OPTIONS]");
    // if (argc >= 3)
    //   app.footer(description);

    // clang-format off
    app.add_option("-i,--input", infile, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "yield output file")
      ->required()
      ->option_text("FILE");
    app.add_option("-s,--step", step_size, "step size in extrapolations")
      ->default_val(step_size);
    app.add_option("-r,--seed", seed, "seed for random number generator")
      ->default_val(seed);
    app.add_flag("-p,--paired-end", paired_end, "input is paired end read file");
    app.add_flag("-v,--verbose", verbose, "print more info");
    // clang-format on

    if (argc < 3) {
      // std::println("{}", app.help());
      std::cout << app.help() << '\n';
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const auto input_format = get_input_format_type(infile);
    if (is_unknown(input_format)) {
      std::cerr << "unknown input format\n";
      return EXIT_FAILURE;
    }

    const auto [n_reads, counts_hist] = [&] {
      switch (input_format) {
      case input_format_type::hist:
        return load_histogram(infile);
      case input_format_type::counts:
        return load_counts(infile);
#ifdef HAVE_HTSLIB
      case input_format_type::bam:
        return paired_end ? load_counts_BAM_pe(n_threads, infile)
                          : load_counts_BAM_se(n_threads, infile);
#endif
      default:  // case input_format_type::vals:
        return paired_end ? load_counts_bed_pe(infile)
                          : load_counts_bed_se(infile);
      }
    }();

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

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open output file: " + outfile);

    out << "total_reads" << '\t' << "distinct_reads\n"
        << 0 << '\t' << 0 << '\n';
    for (std::size_t i = step_size; i <= upper_limit; i += step_size)
      out << i << '\t'
          << interpolate_distinct(counts_hist, total_reads, distinct_reads, i)
          << '\n';
  }
  catch (const std::exception &e) {
    std::cerr << "ERROR:\t" << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
