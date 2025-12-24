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

#include "pop_size.hpp"
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

auto
pop_size::main(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
  try {
    static constexpr auto max_iter_per_bootstrap = 100;
    static constexpr auto default_max_extrap = 1000000000ul;
    static constexpr auto min_required_counts = 4ul;
    static constexpr auto min_required_counts_error_message =
      "max count before zero is less than min required count (4)";

    std::string outfile;
    std::string input_file_name;
    std::string histogram_outfile;

    // NOLINTBEGIN(*-avoid-magic-numbers)
    std::size_t orig_max_terms = 100;
    double max_extrap{};
    std::size_t n_desired_steps = 50;
    std::size_t n_bootstraps = 100;
    int diagonal = 0;
    double c_level = 0.95;
    std::uint32_t seed = 408;
    // NOLINTEND(*-avoid-magic-numbers)

    // flags
    bool verbose{false};
    bool paired_end{false};
    bool single_estimate{false};
    bool allow_defects{false};

#ifdef HAVE_HTSLIB
    std::uint32_t n_threads{1};
#endif

    CLI::App app{rlstrip(pop_size::about_msg)};
    argv = app.ensure_utf8(argv);
    app.usage("\nUsage: preseq pop_size [OPTIONS]");
    if (argc >= 2)
      app.footer(rlstrip(pop_size::footer_msg));

    // clang-format off
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "output file")
      ->option_text("FILE");
    app.add_option("-e,--extrap", max_extrap, "maximum extrapolation");
    app.add_option("-s,--steps", n_desired_steps, "number of steps");
    app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--ci-level", c_level, "level for confidence intervals");
    app.add_option("-x,--terms", orig_max_terms, "maximum terms in estimator");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_flag("-p,--paired-end", paired_end, "input is paired end read file");
    app.add_flag("-q,--quick", single_estimate, "no bootstraps for confidence intervals");
    app.add_flag("-D,--defects", allow_defects, "no testing for defects");
    app.add_flag("-v,--verbose", verbose, "print more info");
    // clang-format on

    if (argc < 3) {
      // std::println("{}", app.help());
      std::cout << app.help() << '\n';
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const auto input_format = get_input_format_type(input_file_name);
    if (is_unknown(input_format)) {
      std::cerr << "unknown input format\n";
      return EXIT_FAILURE;
    }

    if (verbose) {
      std::cerr << "INPUT FORMAT: " << to_string(input_format) << '\n';
      if (is_bam(input_format) || is_bed(input_format))
        std::cerr << "PAIRED END: " << std::boolalpha << paired_end << '\n';
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
    const auto distinct_reads =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::size_t first_zero = 1;
    while (first_zero < std::size(counts_hist) && counts_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = std::min(orig_max_terms, first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

    if (max_extrap == 0.0)
      max_extrap = default_max_extrap * distinct_reads;

    const auto step_size =
      (max_extrap - distinct_reads) / static_cast<double>(n_desired_steps);

    const std::size_t distinct_counts =
      std::count_if(std::cbegin(counts_hist), std::cend(counts_hist),
                    [](const double x) { return x > 0.0; });

    if (verbose)
      std::cerr << "TOTAL READS     = " << n_reads << '\n'
                << "DISTINCT READS  = " << distinct_reads << '\n'
                << "DISTINCT COUNTS = " << distinct_counts << '\n'
                << "MAX COUNT       = " << max_observed_count << '\n'
                << "COUNTS OF 1     = " << counts_hist[1] << '\n'
                << "MAX TERMS       = " << orig_max_terms << '\n';

    if (!histogram_outfile.empty())
      report_histogram(histogram_outfile, counts_hist);

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(counts_hist);
    if (two_fold_extrap < 0.0)
      throw std::runtime_error(
        "Saturation expected at double initial sample size."
        " Unable to extrapolate");

    // check that min required count is satisfied
    if (orig_max_terms < min_required_counts)
      throw std::runtime_error(min_required_counts_error_message);

    if (verbose)
      std::cerr << "[ESTIMATING YIELD CURVE]\n";

    std::vector<double> yield_estimates;

    if (single_estimate) {
      const bool success = extrap_single_estimate(
        verbose, allow_defects, counts_hist, orig_max_terms, diagonal,
        step_size, max_extrap, yield_estimates);
      if (!success)
        throw std::runtime_error(
          "single estimate failed, run full mode for estimates");

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile.data());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS\tEXPECTED_DISTINCT\n";
      out << 0 << '\t' << 0 << '\n';
      for (auto i = 0u; i < std::size(yield_estimates); ++i)
        out << (i + 1) * step_size << '\t' << yield_estimates[i] << '\n';
    }
    else {
      if (verbose)
        std::cerr << "[BOOTSTRAPPING HISTOGRAM]\n";

      const std::size_t max_iter = max_iter_per_bootstrap * n_bootstraps;

      std::vector<std::vector<double>> bootstrap_estimates;
      extrap_bootstrap(verbose, allow_defects, seed, counts_hist, n_bootstraps,
                       orig_max_terms, diagonal, step_size, max_extrap,
                       max_iter, bootstrap_estimates);

      if (verbose)
        std::cerr << "[COMPUTING CONFIDENCE INTERVALS]\n";
      // yield ci
      std::vector<double> yield_upper_ci_lognorm, yield_lower_ci_lognorm;

      vector_median_and_ci(bootstrap_estimates, c_level, yield_estimates,
                           yield_lower_ci_lognorm, yield_upper_ci_lognorm);
      if (verbose)
        std::cerr << "[WRITING OUTPUT]\n";

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile);
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      if (!outfile.empty() && !out)
        throw std::runtime_error("failed to open outfile file: " + outfile);

      const std::size_t n_ests = std::size(yield_estimates) - 1;
      if (n_ests < 2)
        throw std::runtime_error(
          "problem with number of estimates in pop_size");

      const bool converged =
        (yield_estimates[n_ests] - yield_estimates[n_ests - 1] < 1.0);

      out << "pop_size_estimate" << '\t' << "lower_ci" << '\t' << "upper_ci"
          << '\n';
      out << yield_estimates.back() << '\t' << yield_lower_ci_lognorm.back()
          << '\t' << yield_upper_ci_lognorm.back();
      if (!converged)
        out << "\tnot_converged";
      out << '\n';
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
