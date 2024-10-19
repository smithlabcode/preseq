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

#include "bound_pop.hpp"

#include "common.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

#include <OptionParser.hpp>

#include <unistd.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using std::cerr;
using std::endl;
using std::isfinite;
using std::min;
using std::mt19937;
using std::runtime_error;
using std::setprecision;
using std::string;
using std::vector;

namespace fs = std::filesystem;

// BOUND_UNOBS: bounding n_0
int
bound_pop_main(const int argc, const char *argv[]) {
  try {
    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    bool QUICK_MODE = false;

    string outfile;

#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
    uint32_t n_threads{1};
#endif

    size_t max_num_points = 10;
    double tolerance = 1e-20;
    size_t n_bootstraps = 500;
    double c_level = 0.95;
    size_t max_iter = 100;
    uint64_t seed = 408;

    const string description = R"(
Estimate the size of the underlying population based on counts  \
of observed species in an initial sample.)";

    /********** GET COMMAND LINE ARGUMENTS FOR BOUND_POP ***********/
    OptionParser opt_parse(fs::path(argv[1]).filename(), description,
                           "<input-file>");
    opt_parse.add_opt("output", 'o',
                      "species richness output file "
                      "(default: stdout)",
                      false, outfile);
    opt_parse.add_opt("max_num_points", 'p',
                      "maximum number of points in "
                      "quadrature estimates",
                      false, max_num_points);
    opt_parse.add_opt("tolerance", 't', "numerical tolerance", false,
                      tolerance);
    opt_parse.add_opt("bootstraps", 'n', "number of bootstraps", false,
                      n_bootstraps);
    opt_parse.add_opt("clevel", 'c', "level for confidence intervals", false,
                      c_level);
    opt_parse.add_opt("verbose", 'v', "print more information", false, VERBOSE);
    opt_parse.add_opt("pe", 'P', "input is paired end read file", false,
                      PAIRED_END);
    opt_parse.add_opt("hist", 'H',
                      "input is a text file containing the "
                      "observed histogram",
                      false, HIST_INPUT);
    opt_parse.add_opt("vals", 'V',
                      "input is a text file containing only the "
                      "observed duplicate counts",
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
    opt_parse.add_opt("quick", 'Q',
                      "quick mode, estimate without bootstrapping", false,
                      QUICK_MODE);
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
    // ****************************************************************

    vector<double> counts_hist;
    size_t n_obs = 0;

    // LOAD VALUES
    if (HIST_INPUT) {
      if (VERBOSE)
        cerr << "HIST_INPUT" << endl;
      n_obs = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (VERBOSE)
        cerr << "VALS_INPUT" << endl;
      n_obs = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (VERBOSE)
        cerr << "PAIRED_END_BAM_INPUT" << endl;
      n_obs = load_counts_BAM_pe(n_threads, input_file_name, counts_hist);
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        cerr << "BAM_INPUT" << endl;
      n_obs = load_counts_BAM_se(n_threads, input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (VERBOSE)
        cerr << "PAIRED_END_BED_INPUT" << endl;
      n_obs = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_obs = load_counts_BED_se(input_file_name, counts_hist);
    }

    const double distinct_obs =
      accumulate(begin(counts_hist), end(counts_hist), 0.0);

    vector<double> measure_moments;
    // mu_r = (r + 1)! n_{r+1} / n_1
    size_t idx = 1;
    while (counts_hist[idx] > 0 && idx <= counts_hist.size()) {
      // idx + 1 because function calculates (x-1)!
      measure_moments.push_back(
        exp(factorial(idx + 1) + log(counts_hist[idx]) - log(counts_hist[1])));
      if (!isfinite(measure_moments.back())) {
        measure_moments.pop_back();
        break;
      }
      idx++;
    }

    if (VERBOSE) {
      cerr << "TOTAL OBSERVATIONS     = " << n_obs << endl
           << "DISTINCT OBSERVATIONS  = " << distinct_obs << endl
           << "MAX COUNT              = " << counts_hist.size() - 1 << endl;

      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0)
          cerr << i << '\t' << setprecision(16) << counts_hist[i] << endl;

      cerr << "OBSERVED MOMENTS" << endl;
      for (size_t i = 0; i < measure_moments.size(); i++)
        cerr << std::setprecision(16) << measure_moments[i] << endl;
    }

    if (QUICK_MODE) {
      if (measure_moments.size() < 2 * max_num_points)
        max_num_points = static_cast<size_t>(floor(measure_moments.size() / 2));
      else
        measure_moments.resize(2 * max_num_points);
      size_t n_points = 0;
      n_points = ensure_pos_def_mom_seq(measure_moments, tolerance, VERBOSE);
      if (VERBOSE)
        cerr << "n_points = " << n_points << endl;

      MomentSequence obs_mom_seq(measure_moments);

      if (VERBOSE) {
        for (size_t k = 0; k < obs_mom_seq.alpha.size(); k++)
          cerr << "alpha_" << k << '\t';
        cerr << endl;
        for (size_t k = 0; k < obs_mom_seq.alpha.size(); k++)
          cerr << obs_mom_seq.alpha[k] << '\t';
        cerr << endl;

        for (size_t k = 0; k < obs_mom_seq.beta.size(); k++)
          cerr << "beta_" << k << '\t';
        cerr << endl;
        for (size_t k = 0; k < obs_mom_seq.beta.size(); k++)
          cerr << obs_mom_seq.beta[k] << '\t';
        cerr << endl;
      }

      vector<double> points, weights;
      obs_mom_seq.Lower_quadrature_rules(n_points, tolerance, max_iter, points,
                                         weights);

      // renormalize if needed
      const double weights_sum = accumulate(begin(weights), end(weights), 0.0);
      if (weights_sum != 1.0)
        for (size_t i = 0; i < weights.size(); i++)
          weights[i] = weights[i] / weights_sum;

      if (VERBOSE) {
        cerr << "points = " << endl;
        for (size_t i = 0; i < points.size(); i++)
          cerr << points[i] << '\t';
        cerr << endl;

        cerr << "weights = " << endl;
        for (size_t i = 0; i < weights.size(); i++)
          cerr << weights[i] << '\t';
        cerr << endl;
      }

      double estimated_unobs = 0.0;

      for (size_t i = 0; i < weights.size(); i++)
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

      out << "quadrature_estimated_unobs" << '\t' << "n_points" << endl;
      out << estimated_unobs << '\t' << n_points << endl;
    }
    // NOT QUICK MODE, BOOTSTRAP
    else {
      vector<double> quad_estimates;

      // setup rng
      srand(time(0) + getpid());
      mt19937 rng(seed);

      // hist may be sparse, to speed up bootstrapping
      // sample only from positive entries
      vector<size_t> counts_hist_distinct_counts;
      vector<double> distinct_counts_hist;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0) {
          counts_hist_distinct_counts.push_back(i);
          distinct_counts_hist.push_back(counts_hist[i]);
        }

      for (size_t iter = 0;
           iter < max_iter && quad_estimates.size() < n_bootstraps; ++iter) {
        if (VERBOSE)
          cerr << "iter=" << "\t" << iter << endl;

        vector<double> sample_hist;
        resample_hist(rng, counts_hist_distinct_counts, distinct_counts_hist,
                      sample_hist);

        const double sampled_distinct =
          accumulate(begin(sample_hist), end(sample_hist), 0.0);

        // initialize moments, 0th moment is 1
        vector<double> bootstrap_moments(1, 1.0);
        // moments[r] = (r + 1)! n_{r+1} / n_1
        for (size_t i = 0; i < 2 * max_num_points; i++) {
          bootstrap_moments.push_back(exp(
            factorial(i + 3) + log(sample_hist[i + 2]) - log(sample_hist[1])));
        }

        size_t n_points = 0;
        n_points =
          ensure_pos_def_mom_seq(bootstrap_moments, tolerance, VERBOSE);
        n_points = min(n_points, max_num_points);
        if (VERBOSE)
          cerr << "n_points = " << n_points << endl;

        MomentSequence bootstrap_mom_seq(bootstrap_moments);

        vector<double> points, weights;
        bootstrap_mom_seq.Lower_quadrature_rules(n_points, tolerance, max_iter,
                                                 points, weights);

        // renormalize if needed
        const double weights_sum =
          accumulate(begin(weights), end(weights), 0.0);
        if (weights_sum != 1.0)
          for (size_t i = 0; i < weights.size(); i++)
            weights[i] = weights[i] / weights_sum;

        double estimated_unobs = 0.0;

        for (size_t i = 0; i < weights.size(); i++)
          estimated_unobs += counts_hist[1] * weights[i] / points[i];

        if (estimated_unobs > 0.0)
          estimated_unobs += sampled_distinct;
        else {
          estimated_unobs = sampled_distinct;
          n_points = 0;
        }

        if (VERBOSE) {
          cerr << "bootstrapped_moments=" << endl;
          for (size_t i = 0; i < bootstrap_moments.size(); i++)
            cerr << bootstrap_moments[i] << endl;
        }
        if (VERBOSE) {
          for (size_t k = 0; k < bootstrap_mom_seq.alpha.size(); k++)
            cerr << "alpha_" << k << '\t';
          cerr << endl;
          for (size_t k = 0; k < bootstrap_mom_seq.alpha.size(); k++)
            cerr << bootstrap_mom_seq.alpha[k] << '\t';
          cerr << endl;

          for (size_t k = 0; k < bootstrap_mom_seq.beta.size(); k++)
            cerr << "beta_" << k << '\t';
          cerr << endl;
          for (size_t k = 0; k < bootstrap_mom_seq.beta.size(); k++)
            cerr << bootstrap_mom_seq.beta[k] << '\t';
          cerr << endl;
        }
        if (VERBOSE) {
          cerr << "points=" << "\t";
          for (size_t i = 0; i < points.size(); i++)
            cerr << points[i] << "\t";
          cerr << endl;
          cerr << "weights=" << "\t";
          for (size_t i = 0; i < weights.size(); i++)
            cerr << weights[i] << "\t";
          cerr << endl;
          cerr << "estimated_unobs=" << "\t" << estimated_unobs << endl;
        }

        quad_estimates.push_back(estimated_unobs);
      }

      double median_estimate, lower_ci, upper_ci;
      median_and_ci(quad_estimates, c_level, median_estimate, lower_ci,
                    upper_ci);

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile);
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "median_estimated_unobs" << '\t' << "lower_ci" << '\t'
          << "upper_ci" << endl;
      out << median_estimate << '\t' << lower_ci << '\t' << upper_ci << endl;
    }
  }
  catch (runtime_error &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
