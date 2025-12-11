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

const auto about_msg = R"(
preseq bound_pop: Estimate the size of the underlying population based on
counts of observed species in an initial sample.
)";

#include "bound_pop.hpp"

#include "common.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

#include "CLI11/CLI11.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>  // IWYU pragma: keep
#include <numeric>
#include <random>
#include <string>
#include <vector>

// NOLINTBEGIN(*-avoid-magic-numbers,*-narrowing-conversions)

static void
report_bootstrapped_moments(const std::vector<double> &bootstrap_moments,
                            const MomentSequence &bootstrap_mom_seq,
                            const std::vector<double> &points,
                            const std::vector<double> &weights,
                            const double estimated_unobs) {
  std::cerr << "bootstrapped_moments=\n";
  for (const auto m : bootstrap_moments)
    std::cerr << m << '\n';
  for (std::size_t k = 0; k < std::size(bootstrap_mom_seq.alpha); k++)
    std::cerr << "alpha_" << k << '\t';
  std::cerr << '\n';
  for (const auto a : bootstrap_mom_seq.alpha)
    std::cerr << a << '\t';
  std::cerr << '\n';
  for (std::size_t k = 0; k < std::size(bootstrap_mom_seq.beta); k++)
    std::cerr << "beta_" << k << '\t';
  std::cerr << '\n';
  for (const auto b : bootstrap_mom_seq.beta)
    std::cerr << b << '\t';
  std::cerr << '\n';
  std::cerr << "points=\t";
  for (const auto p : points)
    std::cerr << p << '\t';
  std::cerr << '\n';
  std::cerr << "weights=\t";
  for (const auto w : weights)
    std::cerr << w << '\t';
  std::cerr << '\n';
  std::cerr << "estimated_unobs=\t" << estimated_unobs << '\n';
}

// bounding n_0
int
bound_pop_main(int argc, char *argv[]) {  // NOLINT (*-avoid-c-arrays)
  try {
    bool verbose = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    bool QUICK_MODE = false;

    std::string input_file_name;
    std::string outfile;
    std::string histogram_outfile;

#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    std::size_t MAX_SEGMENT_LENGTH = 5000;
    std::uint32_t n_threads{1};
#endif

    std::size_t max_num_points = 10;
    double tolerance = 1e-20;
    std::size_t n_bootstraps = 500;
    double c_level = 0.95;
    std::size_t max_iter = 100;
    std::uint32_t seed = 408;

    CLI::App app{rlstrip(about_msg)};
    argv = app.ensure_utf8(argv);
    app.formatter(std::make_shared<preseq_formatter>());
    app.usage("\nUsage: preseq bound_pop [OPTIONS]");
    // if (argc >= 3)
    //   app.footer(rlstrip(description));

    // clang-format off
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile,
                   "species richness output file (default: print to screen)");
    app.add_option("-p,--max_num_points", max_num_points,
                   "maximum number of points in quadrature estimates");
    app.add_option("-t,--tolerance", tolerance, "numerical tolerance");
    app.add_option("-n,--bootstraps", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--clevel", c_level, "level for confidence intervals");
    app.add_flag("-P,--pe", PAIRED_END, "input is paired end read file");
    app.add_flag("-H,--hist", HIST_INPUT,
                   "input is a text file containing the observed histogram");
    app.add_flag("-V,--vals", VALS_INPUT,
                   "input is a text file containing only the observed duplicate counts");
#ifdef HAVE_HTSLIB
    app.add_flag("-B,--bam", BAM_FORMAT_INPUT,
                   "input is in BAM format");
    app.add_option("-l,--seg_len", MAX_SEGMENT_LENGTH,
                   "maximum segment length when merging paired end bam reads");
#endif
    app.add_flag("-Q,--quick", QUICK_MODE,
                   "quick mode, estimate without bootstrapping");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_flag("-v,--verbose", verbose, "print more info");
    // clang-format on

    if (argc < 3) {
      // std::println("{}", app.help());
      std::cout << app.help() << '\n';
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    std::vector<double> counts_hist;
    std::size_t n_obs = 0;

    // LOAD VALUES
    if (HIST_INPUT) {
      if (verbose)
        std::cerr << "HIST_INPUT\n";
      n_obs = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (verbose)
        std::cerr << "VALS_INPUT\n";
      n_obs = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (verbose)
        std::cerr << "PAIRED_END_BAM_INPUT\n";
      n_obs = load_counts_BAM_pe(n_threads, input_file_name, counts_hist);
    }
    else if (BAM_FORMAT_INPUT) {
      if (verbose)
        std::cerr << "BAM_INPUT\n";
      n_obs = load_counts_BAM_se(n_threads, input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (verbose)
        std::cerr << "PAIRED_END_BED_INPUT\n";
      n_obs = load_counts_bed_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (verbose)
        std::cerr << "BED_INPUT\n";
      n_obs = load_counts_bed_se(input_file_name, counts_hist);
    }

    const double distinct_obs =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    std::vector<double> measure_moments;
    // mu_r = (r + 1)! n_{r+1} / n_1
    std::size_t idx = 1;
    while (idx < std::size(counts_hist) && counts_hist[idx]) {
      // idx + 1 because function calculates (x-1)!
      measure_moments.push_back(std::exp(factorial(idx + 1) +
                                         std::log(counts_hist[idx]) -
                                         std::log(counts_hist[1])));
      if (!std::isfinite(measure_moments.back())) {
        measure_moments.pop_back();
        break;
      }
      ++idx;
    }

    if (verbose) {
      std::cerr << "TOTAL OBSERVATIONS     = " << n_obs << '\n'
                << "DISTINCT OBSERVATIONS  = " << distinct_obs << '\n'
                << "MAX COUNT              = " << std::size(counts_hist) - 1
                << '\n';

      std::cerr << "OBSERVED MOMENTS\n";
      for (std::size_t i = 0; i < std::size(measure_moments); i++)
        std::cerr << std::setprecision(16) << measure_moments[i] << '\n';
    }

    if (!histogram_outfile.empty())
      report_histogram(histogram_outfile, counts_hist);

    if (QUICK_MODE) {
      if (std::size(measure_moments) > 2 * max_num_points)
        measure_moments.resize(2 * max_num_points);

      std::size_t n_points =
        ensure_pos_def_mom_seq(measure_moments, tolerance, verbose);
      if (verbose)
        std::cerr << "n_points = " << n_points << '\n';

      MomentSequence obs_mom_seq(measure_moments);

      if (verbose) {
        for (std::size_t k = 0; k < std::size(obs_mom_seq.alpha); k++)
          std::cerr << "alpha_" << k << '\t';
        std::cerr << '\n';
        for (std::size_t k = 0; k < std::size(obs_mom_seq.alpha); k++)
          std::cerr << obs_mom_seq.alpha[k] << '\t';
        std::cerr << '\n';

        for (std::size_t k = 0; k < std::size(obs_mom_seq.beta); k++)
          std::cerr << "beta_" << k << '\t';
        std::cerr << '\n';
        for (std::size_t k = 0; k < std::size(obs_mom_seq.beta); k++)
          std::cerr << obs_mom_seq.beta[k] << '\t';
        std::cerr << '\n';
      }

      std::vector<double> points, weights;
      obs_mom_seq.Lower_quadrature_rules(n_points, tolerance, max_iter, points,
                                         weights);

      // renormalize if needed
      const double weights_sum =
        std::accumulate(std::cbegin(weights), std::cend(weights), 0.0);
      if (weights_sum != 1.0)
        for (std::size_t i = 0; i < std::size(weights); i++)
          weights[i] = weights[i] / weights_sum;

      if (verbose) {
        std::cerr << "points = \n";
        for (const auto p : points)
          std::cerr << p << '\t';
        std::cerr << '\n';

        std::cerr << "weights = \n";
        for (const auto w : weights)
          std::cerr << w << '\t';
        std::cerr << '\n';
      }

      double estimated_unobs = 0.0;

      for (std::size_t i = 0; i < std::size(weights); i++)
        estimated_unobs += counts_hist[1] * weights[i] / points[i];

      if (estimated_unobs > 0.0)
        estimated_unobs += distinct_obs;
      else {
        estimated_unobs = distinct_obs;
        n_points = 0;
      }

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile);
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "quadrature_estimated_unobs" << '\t' << "n_points\n"
          << estimated_unobs << '\t' << n_points << '\n';
    }
    // NOT QUICK MODE, BOOTSTRAP
    else {
      std::vector<double> quad_estimates;

      // setup rng
      std::mt19937 rng(seed);

      // hist may be sparse, to speed up bootstrapping
      // sample only from positive entries
      std::vector<std::size_t> counts_hist_distinct_counts;
      std::vector<double> distinct_counts_hist;
      for (std::size_t i = 0; i < std::size(counts_hist); i++)
        if (counts_hist[i] > 0) {
          counts_hist_distinct_counts.push_back(i);
          distinct_counts_hist.push_back(counts_hist[i]);
        }

      for (std::size_t iter = 0;
           iter < max_iter && std::size(quad_estimates) < n_bootstraps;
           ++iter) {
        if (verbose)
          std::cerr << "iter=" << "\t" << iter << '\n';

        std::vector<double> sample_hist;
        resample_hist(rng, counts_hist_distinct_counts, distinct_counts_hist,
                      sample_hist);

        const double sampled_distinct = std::accumulate(
          std::cbegin(sample_hist), std::cend(sample_hist), 0.0);

        // initialize moments, 0th moment is 1
        std::vector<double> bootstrap_moments(1, 1.0);
        // moments[r] = (r + 1)! n_{r+1} / n_1
        for (std::size_t i = 0; i < 2 * max_num_points; i++)
          bootstrap_moments.push_back(std::exp(factorial(i + 3) +
                                               std::log(sample_hist[i + 2]) -
                                               std::log(sample_hist[1])));

        std::size_t n_points =
          ensure_pos_def_mom_seq(bootstrap_moments, tolerance, verbose);
        n_points = std::min(n_points, max_num_points);
        if (verbose)
          std::cerr << "n_points = " << n_points << '\n';

        MomentSequence bootstrap_mom_seq(bootstrap_moments);

        std::vector<double> points;
        std::vector<double> weights;
        bootstrap_mom_seq.Lower_quadrature_rules(n_points, tolerance, max_iter,
                                                 points, weights);

        // renormalize if needed
        const double weights_sum =
          std::accumulate(std::cbegin(weights), std::cend(weights), 0.0);
        if (weights_sum != 1.0)
          for (std::size_t i = 0; i < std::size(weights); i++)
            weights[i] = weights[i] / weights_sum;

        double estimated_unobs = 0.0;

        for (std::size_t i = 0; i < std::size(weights); i++)
          estimated_unobs += counts_hist[1] * weights[i] / points[i];

        if (estimated_unobs > 0.0)
          estimated_unobs += sampled_distinct;
        else
          estimated_unobs = sampled_distinct;

        if (verbose)
          report_bootstrapped_moments(bootstrap_moments, bootstrap_mom_seq,
                                      points, weights, estimated_unobs);

        quad_estimates.push_back(estimated_unobs);
      }

      double median_estimate{};
      double lower_ci{};
      double upper_ci{};
      median_and_ci(quad_estimates, c_level, median_estimate, lower_ci,
                    upper_ci);

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile);
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "median_estimated_unobs" << '\t' << "lower_ci" << '\t'
          << "upper_ci\n"
          << median_estimate << '\t' << lower_ci << '\t' << upper_ci << '\n';
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-avoid-magic-numbers,*-narrowing-conversions)
