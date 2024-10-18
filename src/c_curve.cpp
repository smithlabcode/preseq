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

#include "continued_fraction.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

#include <GenomicRegion.hpp>
#include <OptionParser.hpp>
#include <smithlab_os.hpp>
#include <smithlab_utils.hpp>

#include <sys/types.h>
#include <unistd.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_map>
#include <vector>

using std::cerr;
using std::endl;
using std::isfinite;
using std::max;
using std::min;
using std::mt19937;
using std::runtime_error;
using std::setprecision;
using std::string;
using std::to_string;
using std::uint64_t;
using std::unordered_map;
using std::vector;

template <typename T>
T
get_counts_from_hist(const vector<T> &h) {
  T c = 0.0;
  for (size_t i = 0; i < h.size(); ++i)
    c += i * h[i];
  return c;
}

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

template <typename T>
T
quantile_from_sorted_vector(const vector<T> sorted_data, const size_t stride,
                            const size_t n, const double f) {
  const double index = f * (n - 1);
  const size_t lhs = static_cast<int>(index);
  const double delta = index - lhs;

  if (n == 0 || sorted_data.empty())
    return 0.0;

  if (lhs == n - 1)
    return sorted_data[lhs * stride];

  return (1 - delta) * sorted_data[lhs * stride] +
         delta * sorted_data[(lhs + 1) * stride];
}

// Confidence interval stuff
static void
median_and_ci(vector<double> estimates,  // by val so we can sort them
              const double ci_level, double &median_estimate,
              double &lower_ci_estimate, double &upper_ci_estimate) {
  assert(!estimates.empty());

  sort(begin(estimates), end(estimates));

  const double alpha = 1.0 - ci_level;
  const size_t N = estimates.size();

  median_estimate = median_from_sorted_vector(estimates, 1, N);
  lower_ci_estimate = quantile_from_sorted_vector(estimates, 1, N, alpha / 2);
  upper_ci_estimate =
    quantile_from_sorted_vector(estimates, 1, N, 1.0 - alpha / 2);
}

static void
vector_median_and_ci(const vector<vector<double>> &bootstrap_estimates,
                     const double ci_level, vector<double> &yield_estimates,
                     vector<double> &lower_ci_lognorm,
                     vector<double> &upper_ci_lognorm) {
  yield_estimates.clear();
  lower_ci_lognorm.clear();
  upper_ci_lognorm.clear();
  assert(!bootstrap_estimates.empty());

  const size_t n_est = bootstrap_estimates.size();
  vector<double> estimates_row(n_est, 0.0);
  for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {
    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];

    double median_estimate, lower_ci_estimate, upper_ci_estimate;
    median_and_ci(estimates_row, ci_level, median_estimate, lower_ci_estimate,
                  upper_ci_estimate);
    sort(begin(estimates_row), end(estimates_row));

    yield_estimates.push_back(median_estimate);
    lower_ci_lognorm.push_back(lower_ci_estimate);
    upper_ci_lognorm.push_back(upper_ci_estimate);
  }
}

template <typename uint_type>
void
multinomial(mt19937 &gen, const vector<double> &mult_probs, uint_type trials,
            vector<uint_type> &result) {
  typedef std::binomial_distribution<uint32_t> binom_dist;

  result.clear();
  result.resize(mult_probs.size());

  double remaining_prob = accumulate(begin(mult_probs), end(mult_probs), 0.0);

  auto r(begin(result));
  auto p(begin(mult_probs));

  while (p != end(mult_probs)) {  // iterate to sample for each category
    *r = binom_dist(trials, (*p) / remaining_prob)(gen);  // take the sample

    remaining_prob -= *p++;  // update remaining probability mass
    trials -= *r++;          // update remaining trials needed
  }

  if (trials > 0)
    throw runtime_error("multinomial sampling failed");
}

// Lanczos approximation for gamma function for x >= 0.5 - essentially an
// approximation for (x-1)!
static double
factorial(double x) {
  // constants
  double LogRootTwoPi = 0.9189385332046727;
  double Euler = 2.71828182845904523536028747135;

  // Approximation for factorial is actually x-1
  x -= 1.0;

  vector<double> lanczos{0.99999999999980993227684700473478,
                         676.520368121885098567009190444019,
                         -1259.13921672240287047156078755283,
                         771.3234287776530788486528258894,
                         -176.61502916214059906584551354,
                         12.507343278686904814458936853,
                         -0.13857109526572011689554707,
                         9.984369578019570859563e-6,
                         1.50563273514931155834e-7};

  double Ag = lanczos[0];

  for (size_t k = 1; k < lanczos.size(); k++)
    Ag += lanczos[k] / (x + k);

  double term1 = (x + 0.5) * log((x + 7.5) / Euler);
  double term2 = LogRootTwoPi + log(Ag);

  return term1 + (term2 - 7.0);
}

// interpolate by explicit calculating the expectation
// for sampling without replacement;
// see K.L Heck 1975
// N total sample size; S the total number of distincts
// n sub sample size
static double
interpolate_distinct(const vector<double> &hist, const size_t N, const size_t S,
                     const size_t n) {
  const double denom =
    factorial(N + 1) - factorial(n + 1) - factorial(N - n + 1);

  vector<double> numer(hist.size(), 0);
  for (size_t i = 1; i < hist.size(); i++) {
    // N - i -n + 1 should be greater than 0
    if (N < i + n) {
      numer[i] = 0;
    }
    else {
      const double x =
        (factorial(N - i + 1) - factorial(n + 1) - factorial(N - i - n + 1));
      numer[i] = exp(x - denom) * hist[i];
    }
  }
  return S - accumulate(begin(numer), end(numer), 0);
}

static double
GoodToulmin2xExtrap(const vector<double> &counts_hist) {
  double two_fold_extrap = 0.0;
  for (size_t i = 0; i < counts_hist.size(); i++)
    two_fold_extrap += pow(-1.0, i + 1) * counts_hist[i];
  return two_fold_extrap;
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
      const size_t MAX_READS_TO_HOLD = 5000000;
      size_t n_paired = 0;
      size_t n_mates = 0;
      n_reads =
        load_counts_BAM_pe(VERBOSE, input_file_name, MAX_SEGMENT_LENGTH,
                           MAX_READS_TO_HOLD, n_paired, n_mates, counts_hist);
      if (VERBOSE)
        cerr << "MERGED PAIRED END READS = " << n_paired << endl
             << "MATES PROCESSED = " << n_mates << endl;
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
