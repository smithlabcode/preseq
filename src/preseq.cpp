/*
 *    preseq: to predict properties of genomic sequencing libraries
 *
 *    Copyright (C) 2013-2015 University of Southern California and
 *             Andrew D. Smith and Timothy Daley
 *
 *    Authors: Timothy Daley, Chao Deng, Victoria Helus, and Andrew Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <numeric>
#include <vector>
#include <iomanip>
#include <sys/types.h>
#include <unistd.h>

#include <cstring>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <smithlab_os.hpp>

#include "continued_fraction.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

using std::string;
using std::min;
using std::vector;
using std::endl;
using std::cerr;
using std::isfinite;

using std::setprecision;
using std::unordered_map;
using std::runtime_error;
using std::to_string;
using std::mt19937;
using std::max;

static const string preseq_version = "3.0.2";

template <typename T> T
get_counts_from_hist(const vector<T> &h) {
  T c = 0.0;
  for (size_t i = 0; i < h.size(); ++i)
    c += i*h[i];
  return c;
}

template<typename T> T
median_from_sorted_vector (const vector<T> sorted_data,
                           const size_t stride, const size_t n) {

  if (n == 0 || sorted_data.empty()) return 0.0;

  const size_t lhs = (n - 1) / 2;
  const size_t rhs = n / 2;

  if (lhs == rhs) return sorted_data[lhs * stride];

  return (sorted_data[lhs * stride] + sorted_data[rhs * stride]) / 2.0;
}

template<typename T> T
quantile_from_sorted_vector (const vector<T> sorted_data,
                             const size_t stride, const size_t n,
                             const double f) {
  const double index = f * (n - 1);
  const size_t lhs = (int)index;
  const double delta = index - lhs;

  if (n == 0 || sorted_data.empty()) return 0.0;

  if (lhs == n - 1) return sorted_data[lhs * stride];

  return (1 - delta) * sorted_data[lhs * stride]
    + delta * sorted_data[(lhs + 1) * stride];
}

// Confidence interval stuff
static void
median_and_ci(vector<double> estimates, // by val so we can sort them
              const double ci_level, double &median_estimate,
              double &lower_ci_estimate, double &upper_ci_estimate) {
  assert(!estimates.empty());

  sort(begin(estimates), end(estimates));

  const double alpha = 1.0 - ci_level;
  const size_t N = estimates.size();

  median_estimate = median_from_sorted_vector(estimates, 1, N);
  lower_ci_estimate = quantile_from_sorted_vector(estimates, 1, N, alpha/2);
  upper_ci_estimate = quantile_from_sorted_vector(estimates, 1, N, 1.0 - alpha/2);
}

static void
vector_median_and_ci(const vector<vector<double> > &bootstrap_estimates,
                     const double ci_level,
                     vector<double> &yield_estimates,
                     vector<double> &lower_ci_lognorm,
                     vector<double> &upper_ci_lognorm) {

  yield_estimates.clear();
  lower_ci_lognorm.clear();
  upper_ci_lognorm.clear();
  assert(!bootstrap_estimates.empty());

  const size_t n_est = bootstrap_estimates.size();
  vector<double> estimates_row(n_est, 0.0);
  double prev_estimate = 0;
  for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {

    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];

    double median_estimate, lower_ci_estimate, upper_ci_estimate;
    median_and_ci(estimates_row, ci_level, median_estimate,
                  lower_ci_estimate, upper_ci_estimate);
    sort(begin(estimates_row), end(estimates_row));

    if(median_estimate - prev_estimate < 1.0)
      break;

    yield_estimates.push_back(median_estimate);
    lower_ci_lognorm.push_back(lower_ci_estimate);
    upper_ci_lognorm.push_back(upper_ci_estimate);

    prev_estimate = median_estimate;
  }
}

template <typename uint_type>
void
multinomial(mt19937 &gen, const vector<double> &mult_probs,
            uint_type trials, vector<uint_type> &result) {

  typedef std::binomial_distribution<unsigned int> binom_dist;

  result.clear();
  result.resize(mult_probs.size());

  double remaining_prob = accumulate(begin(mult_probs), end(mult_probs), 0.0);

  auto r(begin(result));
  auto p(begin(mult_probs));

  while (p != end(mult_probs)) { // iterate to sample for each category

    *r = binom_dist(trials, (*p)/remaining_prob)(gen); // take the sample

    remaining_prob -= *p++; // update remaining probability mass
    trials -= *r++;         // update remaining trials needed
  }

  if (trials > 0)
    throw runtime_error("multinomial sampling failed");
}

// Lanczos approximation for gamma function for x >= 0.5 - essentially an
// approximation for (x-1)!
static double
factorial (double x) {

  // constants
  double LogRootTwoPi = 0.9189385332046727;
  double Euler = 2.71828182845904523536028747135;

  // Approximation for factorial is actually x-1
  x -= 1.0;

  vector<double> lanczos {
                          0.99999999999980993227684700473478,
                          676.520368121885098567009190444019,
                          -1259.13921672240287047156078755283,
                          771.3234287776530788486528258894,
                          -176.61502916214059906584551354,
                          12.507343278686904814458936853,
                          -0.13857109526572011689554707,
                          9.984369578019570859563e-6,
                          1.50563273514931155834e-7
  };

  double Ag = lanczos[0];

  for (size_t k=1; k < lanczos.size(); k++)
    Ag += lanczos[k] / (x + k);

  double term1 = (x + 0.5) * log((x + 7.5) / Euler);
  double term2 = LogRootTwoPi + log(Ag);

  return term1 + (term2 - 7.0);
}

////////////////////////////////////////////////////////////////////////
/////  EXTRAP MODE BELOW HERE

// vals_hist[j] = n_{j} = # (counts = j)
// vals_hist_distinct_counts[k] = kth index j s.t. vals_hist[j] > 0
// stores kth index of vals_hist that is positive
// distinct_counts_hist[k] = vals_hist[vals_hist_distinct_counts[k]]
// stores the kth positive value of vals_hist
void
resample_hist(mt19937 &gen, const vector<size_t> &vals_hist_distinct_counts,
              const vector<double> &distinct_counts_hist,
              vector<double> &out_hist) {

  const size_t hist_size = distinct_counts_hist.size();
  vector<unsigned int> sample_distinct_counts_hist(hist_size, 0);

  unsigned int distinct =
    accumulate(begin(distinct_counts_hist), end(distinct_counts_hist), 0.0);

  multinomial(gen, distinct_counts_hist, distinct,
              sample_distinct_counts_hist);

  out_hist.clear();
  out_hist.resize(vals_hist_distinct_counts.back() + 1, 0.0);
  for (size_t i = 0; i < hist_size; i++)
    out_hist[vals_hist_distinct_counts[i]] = sample_distinct_counts_hist[i];
}

// interpolate by explicit calculating the expectation
// for sampling without replacement;
// see K.L Heck 1975
// N total sample size; S the total number of distincts
// n sub sample size
static double
interpolate_distinct(const vector<double> &hist, const size_t N,
                     const size_t S, const size_t n) {

  const double denom = factorial(N + 1) - factorial(n + 1) - factorial(N - n + 1);

  vector<double> numer(hist.size(), 0);
  for (size_t i = 1; i < hist.size(); i++) {
    // N - i -n + 1 should be greater than 0
    if (N < i + n) {
      numer[i] = 0;
    }
    else {

      const double x = (factorial(N - i + 1) - factorial(n + 1) -
                        factorial(N - i - n + 1));
      numer[i] = exp(x - denom)*hist[i];
    }
  }
  return S - accumulate(begin(numer), end(numer), 0);
}


static void
extrapolate_curve(const ContinuedFraction &the_cf,
                  const double initial_distinct,
                  const double vals_sum,
                  const double initial_sample_size,
                  const double step_size,
                  const double max_sample_size,
                  vector<double> &estimates) {

  double curr_samp_sz = initial_sample_size;
  while (curr_samp_sz < max_sample_size) {
    const double fold = (curr_samp_sz - vals_sum)/vals_sum;
    assert(fold >= 0.0);
    estimates.push_back(initial_distinct + fold*the_cf(fold));
    curr_samp_sz += step_size;
  }
}

void
extrap_bootstrap(const bool VERBOSE, const bool allow_defects,
                 const unsigned long int seed,
                 const vector<double> &orig_hist,
                 const size_t n_bootstraps, const size_t orig_max_terms,
                 const int diagonal, const double bin_step_size,
                 const double max_extrap, const size_t max_iter,
                 vector<vector<double> > &bootstrap_estimates) {
  // clear returning vectors
  bootstrap_estimates.clear();

  //setup rng
  srand(time(0) + getpid());
  mt19937 rng(seed);

  // const double vals_sum = get_counts_from_hist(orig_hist);
  const double initial_distinct =
    accumulate(begin(orig_hist), end(orig_hist), 0.0);

  vector<size_t> orig_hist_distinct_counts;
  vector<double> distinct_orig_hist;
  for (size_t i = 0; i < orig_hist.size(); i++)
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }

  for (size_t iter = 0; (iter < max_iter && bootstrap_estimates.size() < n_bootstraps); ++iter) {

    vector<double> yield_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);

    const double sample_vals_sum = get_counts_from_hist(hist);

    // resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    // compute complexity curve by random sampling w/out replacement
    const size_t distinct = accumulate(begin(hist), end(hist), 0.0);
    size_t curr_sample_sz = bin_step_size;
    while (curr_sample_sz < sample_vals_sum) {
      yield_vector.push_back(interpolate_distinct(hist, sample_vals_sum,
                                                  distinct, curr_sample_sz));
      curr_sample_sz += bin_step_size;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t first_zero = 1;
    while (first_zero < hist.size() && hist[first_zero] > 0)
      ++first_zero;

    size_t max_terms = min(orig_max_terms, first_zero - 1);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 1);

    bool successful_bootstrap = false;
    // defect mode, simple extrapolation
    if (allow_defects) {

      vector<double> ps_coeffs;
      for (size_t j = 1; j <= max_terms; j++)
        ps_coeffs.push_back(hist[j]*std::pow(-1.0, j + 1));

      const ContinuedFraction defect_cf(ps_coeffs, diagonal, max_terms);

      extrapolate_curve(defect_cf, initial_distinct, sample_vals_sum,
                        curr_sample_sz, bin_step_size,
                        max_extrap, yield_vector);
      // no checking of curve in defect mode
      bootstrap_estimates.push_back(yield_vector);
      successful_bootstrap = true;
    }
    else {

      // refit curve for lower bound
      const ContinuedFractionApproximation lower_cfa(diagonal, max_terms);
      const ContinuedFraction lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));

      // extrapolate the curve start
      if (lower_cf.is_valid()) {

        extrapolate_curve(lower_cf, initial_distinct, sample_vals_sum,
                          curr_sample_sz, bin_step_size,
                          max_extrap, yield_vector);
        // sanity check
        if (check_yield_estimates_stability(yield_vector)) {
          bootstrap_estimates.push_back(yield_vector);
          successful_bootstrap = true;
        }
      }
    }
    if (VERBOSE)
      cerr << (successful_bootstrap ? '.' : '_');
  }
  if (VERBOSE)
    cerr << endl;
  if (bootstrap_estimates.size() < n_bootstraps)
    throw runtime_error("too many defects in the approximation, "
                        "consider running in defect mode");
}

static bool
extrap_single_estimate(const bool VERBOSE, const bool allow_defects,
                       vector<double> &hist,
                       size_t max_terms, const int diagonal,
                       const double step_size,
                       const double max_extrap,
                       vector<double> &yield_estimate) {

  yield_estimate.clear();

  const double vals_sum = get_counts_from_hist(hist);
  const double initial_distinct = accumulate(begin(hist), end(hist), 0.0);

  // interpolate complexity curve by random sampling w/out replacement
  size_t upper_limit = static_cast<size_t>(vals_sum);
  size_t step = static_cast<size_t>(step_size);
  size_t sample = static_cast<size_t>(step_size);
  for (; sample < upper_limit; sample += step)
    yield_estimate.push_back(interpolate_distinct(hist, upper_limit,
                                                  initial_distinct, sample));

  // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
  size_t first_zero = 1;
  while (first_zero < hist.size() && hist[first_zero] > 0)
    ++first_zero;

  // Ensure we are not using a zero term
  max_terms = min(max_terms, first_zero - 1);

  // refit curve for lower bound (degree of approx is 1 less than
  // max_terms)
  max_terms = max_terms - (max_terms % 2 == 1);

  if (allow_defects) {

    vector<double> ps_coeffs;
    for (size_t j = 1; j <= max_terms; j++)
      ps_coeffs.push_back(hist[j]*std::pow(-1.0, j + 1));

    const ContinuedFraction defect_cf(ps_coeffs, diagonal, max_terms);

    extrapolate_curve(defect_cf, initial_distinct, vals_sum,
                      sample, step_size, max_extrap, yield_estimate);

    if (VERBOSE)
      cerr << defect_cf << endl;
    // NO FAIL! defect mode doesn't care about failure
  }
  else {

    const ContinuedFractionApproximation lower_cfa(diagonal, max_terms);
    const ContinuedFraction lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));

    // extrapolate curve
    if (lower_cf.is_valid()) {
      extrapolate_curve(lower_cf, initial_distinct, vals_sum,
                        sample, step_size, max_extrap, yield_estimate);
    }
    else {
      // FAIL!
      // lower_cf unacceptable, need to bootstrap to obtain estimates
      return false;
    }

    if (VERBOSE)
      cerr << lower_cf << endl;
  }
  // SUCCESS!!
  return true;
}


static double
GoodToulmin2xExtrap(const vector<double> &counts_hist) {
  double two_fold_extrap = 0.0;
  for (size_t i = 0; i < counts_hist.size(); i++)
    two_fold_extrap += pow(-1.0, i + 1)*counts_hist[i];
  return two_fold_extrap;
}


static void
write_predicted_complexity_curve(const string &outfile,
                                 const double c_level, const double step_size,
                                 const vector<double> &yield_estimates,
                                 const vector<double> &yield_lower_ci_lognorm,
                                 const vector<double> &yield_upper_ci_lognorm) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOWER_" << c_level << "CI\t"
      << "UPPER_" << c_level << "CI" << endl;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < yield_estimates.size(); ++i)
    out << (i + 1)*step_size << '\t'
        << yield_estimates[i] << '\t'
        << yield_lower_ci_lognorm[i] << '\t'
        << yield_upper_ci_lognorm[i] << endl;
}


// ADS: functions same, header different (above and this one)
static void
write_predicted_coverage_curve(const string &outfile,
                               const double c_level,
                               const double base_step_size,
                               const size_t bin_size,
                               const vector<double> &cvrg_estimates,
                               const vector<double> &cvrg_lower_ci_lognorm,
                               const vector<double> &cvrg_upper_ci_lognorm) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "TOTAL_BASES\tEXPECTED_COVERED_BASES\t"
      << "LOWER_" << 100*c_level << "%CI\t"
      << "UPPER_" << 100*c_level << "%CI" << endl;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < cvrg_estimates.size(); ++i)
    out << (i + 1)*base_step_size << '\t'
        << cvrg_estimates[i]*bin_size << '\t'
        << cvrg_lower_ci_lognorm[i]*bin_size << '\t'
        << cvrg_upper_ci_lognorm[i]*bin_size << endl;
}


static int
lc_extrap(const bool pop_size, const int argc, const char **argv) {

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
    unsigned long int seed = 408;

    if (pop_size) {
      // ADS: extrapolate far, without too many steps...
      max_extrap = 1.0e20;
      step_size = 1.0e18;
    }

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

    string description;
    if (!pop_size) {
      description =
        "Extrapolate the complexity of a library. This is the approach   \
        described in Daley & Smith (2013). The method applies rational   \
        function approximation via continued fractions with the          \
        original goal of estimating the number of distinct reads that a  \
        sequencing library would yield upon deeper sequencing. This      \
        method has been used for many different purposes since then.";
    }
    else {
      description =
        "Determine the estimate of the number of unique classes in a library.";
    }

    /********** GET COMMAND LINE ARGUMENTS  FOR LC EXTRAP ***********/

    OptionParser opt_parse(strip_path(argv[1]), description, "<input-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation", false, max_extrap);
    opt_parse.add_opt("step",'s',"extrapolation step size", false, step_size);
    opt_parse.add_opt("boots",'n',"number of bootstraps", false, n_bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals", false, c_level);
    opt_parse.add_opt("terms",'x',"maximum terms in estimator", false, orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more info", false, VERBOSE);
#ifdef HAVE_HTSLIB
    opt_parse.add_opt("bam", 'B', "input is in BAM format",
                      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads",
                      false, MAX_SEGMENT_LENGTH);
#endif
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
                      false, PAIRED_END);
    opt_parse.add_opt("vals", 'V',
                      "input is a text file containing only the observed counts",
                      false, VALS_INPUT);
    opt_parse.add_opt("hist", 'H',
                      "input is a text file containing the observed histogram",
                      false, HIST_INPUT);
    opt_parse.add_opt("quick", 'Q',
                      "quick mode (no bootstraps) for confidence intervals",
                      false, SINGLE_ESTIMATE);
    opt_parse.add_opt("defects", 'D', "no testing for defects", false, allow_defects);
    opt_parse.add_opt("seed", 'r', "seed for random number generator",
                      false, seed);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    // ADS: suspect bug below; "-about" isn't working.
    opt_parse.parse(argc-1, argv+1, leftover_args);
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
      n_reads = load_counts_BAM_pe(VERBOSE, input_file_name,
                                   MAX_SEGMENT_LENGTH,
                                   MAX_READS_TO_HOLD, n_paired,
                                   n_mates, counts_hist);
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
    else { // default is single end bed file
      if (VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }
    /************ done loading input **********************************/

    const size_t max_observed_count = counts_hist.size() - 1;
    const double distinct_reads =
      accumulate(begin(counts_hist), end(counts_hist), 0.0);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t first_zero = 1;
    while (first_zero < counts_hist.size() && counts_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = min(orig_max_terms, first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

    const size_t distinct_counts =
      std::count_if(begin(counts_hist), end(counts_hist),
                    [](const double x) {return x > 0.0;});

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
      throw runtime_error("Saturation expected at double initial sample size."
                          " Unable to extrapolate");

    // const size_t total_reads = get_counts_from_hist(counts_hist);

    //assert(total_reads == n_reads); // ADS: why commented out?

    // check that min required count is satisfied
    if (orig_max_terms < min_required_counts)
      throw runtime_error(min_required_counts_error_message);

    if (VERBOSE)
      cerr << "[ESTIMATING YIELD CURVE]" << endl;
    vector<double> yield_estimates;

    if (SINGLE_ESTIMATE) {

      const bool single_estimate_success =
        extrap_single_estimate(VERBOSE, allow_defects, counts_hist, orig_max_terms,
                               diagonal, step_size, max_extrap, yield_estimates);
      // IF FAILURE, EXIT
      if (!single_estimate_success)
        throw runtime_error("single estimate failed, run "
                            "full mode for estimates");

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS\tEXPECTED_DISTINCT" << endl;
      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << endl;
      for (size_t i = 0; i < yield_estimates.size(); ++i)
        out << (i + 1)*step_size << '\t' << yield_estimates[i] << endl;
    }
    else {
      if (VERBOSE)
        cerr << "[BOOTSTRAPPING HISTOGRAM]" << endl;

      const size_t max_iter = 100*n_bootstraps;

      vector<vector <double> > bootstrap_estimates;
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

      if (!pop_size) {
        write_predicted_complexity_curve(outfile, c_level, step_size,
                                         yield_estimates, yield_lower_ci_lognorm,
                                         yield_upper_ci_lognorm);
      }
      else {
        std::ofstream of;
        if (!outfile.empty()) of.open(outfile.c_str());
        std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

        out.setf(std::ios_base::fixed, std::ios_base::floatfield);
        out.precision(1);

        out << "pop_size_estimate" << '\t'
            << "lower_ci" << '\t' << "upper_ci" << endl;
        out << yield_estimates[yield_estimates.size() - 2] << '\t'
            << yield_lower_ci_lognorm[yield_estimates.size() - 2] << '\t'
            << yield_upper_ci_lognorm[yield_estimates.size() - 2] << endl;
        out << yield_estimates.back() << '\t'
            << yield_lower_ci_lognorm.back() << '\t'
            << yield_upper_ci_lognorm.back() << endl;
      }
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


///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///// GC_EXTRAP: predicting genomic coverage
/////

static int
gc_extrap(const int argc, const char **argv) {

  try {

    const size_t MIN_REQUIRED_COUNTS = 4;

    int diagonal = 0;
    size_t orig_max_terms = 100;
    size_t bin_size = 10;
    bool VERBOSE = false;
    string outfile;
    double base_step_size = 1.0e8;
    size_t max_width = 10000;
    bool SINGLE_ESTIMATE = false;
    double max_extrap = 1.0e12;
    size_t n_bootstraps = 100;
    unsigned long int seed = 408;
    bool allow_defects = false;

    bool NO_SEQUENCE = false;
    double c_level = 0.95;

    const string description =
      "Extrapolate the size of the covered genome by mapped reads. This \
      approach is described in Daley & Smith (2014). The method is the  \
      same as for lc_extrap: using rational function approximation to   \
      a power-series expansion for the number of \"unobserved\" bases   \
      in the initial sample. The gc_extrap method is adapted to deal    \
      with individual nucleotides rather than distinct reads.";

    // ********* GET COMMAND LINE ARGUMENTS  FOR GC EXTRAP **********
    OptionParser opt_parse(strip_path(argv[1]), description,
                           "<input-file>");
    opt_parse.add_opt("output", 'o', "coverage yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("max_width", 'w', "max fragment length, "
                      "set equal to read length for single end reads",
                      false, max_width);
    opt_parse.add_opt("bin_size", 'b', "bin size",
                      false, bin_size);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation in base pairs",
                      false, max_extrap);
    opt_parse.add_opt("step",'s',"step size in bases between extrapolations",
                      false, base_step_size);
    opt_parse.add_opt("bootstraps",'n',"number of bootstraps",
                      false, n_bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals",
                      false, c_level);
    opt_parse.add_opt("terms",'x',"maximum number of terms",
                      false, orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.add_opt("bed", 'B',
                      "input is in bed format without sequence information",
                      false, NO_SEQUENCE);
    opt_parse.add_opt("quick",'Q',
                      "quick mode: run gc_extrap without "
                      "bootstrapping for confidence intervals",
                      false, SINGLE_ESTIMATE);
    opt_parse.add_opt("defects", 'D',
                      "defects mode to extrapolate without testing for defects",
                      false, allow_defects);
    opt_parse.add_opt("seed", 'r', "seed for random number generator",
                      false, seed);
    opt_parse.set_show_defaults();

    vector<string> leftover_args;
    opt_parse.parse(argc-1, argv+1, leftover_args);
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

    vector<double> coverage_hist;
    size_t n_reads = 0;
    if (VERBOSE)
      cerr << "LOADING READS" << endl;

    if (NO_SEQUENCE) {
      if (VERBOSE)
        cerr << "BED FORMAT" << endl;
      n_reads = load_coverage_counts_GR(input_file_name, seed, bin_size,
                                        max_width, coverage_hist);
    }
    else {
      if (VERBOSE)
        cerr << "MAPPED READ FORMAT" << endl;
      n_reads = load_coverage_counts_MR(VERBOSE, input_file_name, seed, bin_size,
                                        max_width, coverage_hist);
    }

    const double total_bins = get_counts_from_hist(coverage_hist);

    const double distinct_bins =
      accumulate(coverage_hist.begin(), coverage_hist.end(), 0.0);

    const double avg_bins_per_read = total_bins/n_reads;
    double bin_step_size = base_step_size/bin_size;

    const size_t max_observed_count = coverage_hist.size() - 1;

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t first_zero = 1;
    while (first_zero < coverage_hist.size() && coverage_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = min(orig_max_terms, first_zero - 1);

    if (VERBOSE)
      cerr << "TOTAL READS         = " << n_reads << endl
           << "BASE STEP SIZE      = " << base_step_size << endl
           << "BIN STEP SIZE       = " << bin_step_size << endl
           << "TOTAL BINS          = " << total_bins << endl
           << "BINS PER READ       = " << avg_bins_per_read << endl
           << "DISTINCT BINS       = " << distinct_bins << endl
           << "TOTAL BASES         = " << total_bins*bin_size << endl
           << "TOTAL COVERED BASES = " << distinct_bins*bin_size << endl
           << "MAX COVERAGE COUNT  = " << max_observed_count << endl
           << "COUNTS OF 1         = " << coverage_hist[1] << endl;

    if (VERBOSE) {
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

    if (VERBOSE)
      cerr << "[ESTIMATING COVERAGE CURVE]" << endl;

    vector<double> coverage_estimates;

    if (SINGLE_ESTIMATE) {

      bool SINGLE_ESTIMATE_SUCCESS =
        extrap_single_estimate(VERBOSE, allow_defects, coverage_hist, orig_max_terms, diagonal,
                               bin_step_size, max_extrap/bin_size,
                               coverage_estimates);
      // IF FAILURE, EXIT
      if (!SINGLE_ESTIMATE_SUCCESS)
        throw runtime_error("SINGLE ESTIMATE FAILED, NEED TO RUN IN "
                            "FULL MODE FOR ESTIMATES");

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_BASES\tEXPECTED_DISTINCT" << endl;

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << endl;
      for (size_t i = 0; i < coverage_estimates.size(); ++i)
        out << (i + 1)*base_step_size << '\t'
            << coverage_estimates[i]*bin_size << endl;
    }
    else {

      if (VERBOSE)
        cerr << "[BOOTSTRAPPING HISTOGRAM]" << endl;

      const size_t max_iter = 10*n_bootstraps;

      vector<vector<double> > bootstrap_estimates;
      extrap_bootstrap(VERBOSE, allow_defects, seed, coverage_hist, n_bootstraps, orig_max_terms,
                       diagonal, bin_step_size, max_extrap/bin_size,
                       max_iter, bootstrap_estimates);

      if (VERBOSE)
        cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      vector<double> coverage_upper_ci_lognorm, coverage_lower_ci_lognorm;
      vector_median_and_ci(bootstrap_estimates, c_level,
                           coverage_estimates, coverage_lower_ci_lognorm,
                           coverage_upper_ci_lognorm);

      if (VERBOSE)
        cerr << "[WRITING OUTPUT]" << endl;

      write_predicted_coverage_curve(outfile, c_level, base_step_size,
                                     bin_size, coverage_estimates,
                                     coverage_lower_ci_lognorm,
                                     coverage_upper_ci_lognorm);

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



////////////////////////////////////////////////////////////////////////
/////  C_CURVE BELOW HERE

static int
c_curve(const int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    unsigned long int seed = 408;

    string outfile;

    size_t upper_limit = 0;
    double step_size = 1e6;

#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    const string description =
      "Generate the complexity curve for data. This does not extrapolate, \
      but instead resamples from the given data.";

    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]), description, "<input-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("step",'s',"step size in extrapolations",
                      false, step_size);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
                      false, PAIRED_END);
    opt_parse.add_opt("hist", 'H',
                      "input is a text file containing the observed histogram",
                      false, HIST_INPUT);
    opt_parse.add_opt("vals", 'V',
                      "input is a text file containing only the observed counts",
                      false, VALS_INPUT);
#ifdef HAVE_HTSLIB
    opt_parse.add_opt("bam", 'B', "input is in BAM format",
                      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads",
                      false, MAX_SEGMENT_LENGTH);
#endif
    opt_parse.add_opt("seed", 'r', "seed for random number generator",
                      false, seed);
    opt_parse.set_show_defaults();

    vector<string> leftover_args;
    opt_parse.parse(argc-1, argv+1, leftover_args);
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
    srand(time(0) + getpid()); //give the random fxn a new seed
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
      n_reads = load_counts_BAM_pe(VERBOSE, input_file_name,
                                   MAX_SEGMENT_LENGTH, MAX_READS_TO_HOLD,
                                   n_paired, n_mates, counts_hist);
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
    else { // default is single end bed file
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
                    [](const double x) {return x > 0.0;});

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
      upper_limit = n_reads; //set upper limit to equal the number of molecules

    //handles output of c_curve
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    //prints the complexity curve
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


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// BOUND_UNOBS: bounding n_0

static int
bound_pop(const int argc, const char **argv) {

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
#endif

    size_t max_num_points = 10;
    double tolerance = 1e-20;
    size_t n_bootstraps = 500;
    double c_level = 0.95;
    size_t max_iter = 100;
    unsigned long int seed = 408;

    const string description =
      "Estimate the size of the underlying population based on counts \
      of observed species in an initial sample.";

    /********** GET COMMAND LINE ARGUMENTS FOR BOUND_POP ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           description, "<input-file>");
    opt_parse.add_opt("output", 'o', "species richness output file "
                      "(default: stdout)", false , outfile);
    opt_parse.add_opt("max_num_points",'p',"maximum number of points in "
                      "quadrature estimates", false, max_num_points);
    opt_parse.add_opt("tolerance", 't', "numerical tolerance",
                      false, tolerance);
    opt_parse.add_opt("bootstraps", 'n', "number of bootstraps",
                      false, n_bootstraps);
    opt_parse.add_opt("clevel", 'c', "level for confidence intervals",
                      false, c_level);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
                      false, PAIRED_END);
    opt_parse.add_opt("hist", 'H', "input is a text file containing the "
                      "observed histogram", false, HIST_INPUT);
    opt_parse.add_opt("vals", 'V', "input is a text file containing only the "
                      "observed duplicate counts", false, VALS_INPUT);
#ifdef HAVE_HTSLIB
    opt_parse.add_opt("bam", 'B', "input is in BAM format",
                      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads",
                      false, MAX_SEGMENT_LENGTH);
#endif
    opt_parse.add_opt("quick", 'Q', "quick mode, estimate without bootstrapping",
                      false, QUICK_MODE);
    opt_parse.add_opt("seed", 'r', "seed for random number generator",
                      false, seed);
    opt_parse.set_show_defaults();

    vector<string> leftover_args;
    opt_parse.parse(argc-1, argv+1, leftover_args);
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
      const size_t MAX_READS_TO_HOLD = 5000000;
      size_t n_paired = 0;
      size_t n_mates = 0;
      n_obs = load_counts_BAM_pe(VERBOSE, input_file_name,
                                 MAX_SEGMENT_LENGTH,
                                 MAX_READS_TO_HOLD, n_paired,
                                 n_mates, counts_hist);
      if (VERBOSE) {
        cerr << "MERGED PAIRED END READS = " << n_paired << endl;
        cerr << "MATES PROCESSED = " << n_mates << endl;
      }
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        cerr << "BAM_INPUT" << endl;
      n_obs = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (VERBOSE)
        cerr << "PAIRED_END_BED_INPUT" << endl;
      n_obs = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else { // default is single end bed file
      if (VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_obs = load_counts_BED_se(input_file_name, counts_hist);
    }

    const double distinct_obs = accumulate(begin(counts_hist), end(counts_hist), 0.0);

    vector<double> measure_moments;
    // mu_r = (r + 1)! n_{r+1} / n_1
    size_t idx = 1;
    while (counts_hist[idx] > 0 && idx <= counts_hist.size()) {
      // idx + 1 because function calculates (x-1)!
      measure_moments.push_back(exp(factorial(idx + 1) +
                                    log(counts_hist[idx]) -
                                    log(counts_hist[1])));
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
      if (measure_moments.size() < 2*max_num_points)
        max_num_points = static_cast<size_t>(floor(measure_moments.size()/2));
      else
        measure_moments.resize(2*max_num_points);
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
      obs_mom_seq.Lower_quadrature_rules(VERBOSE, n_points, tolerance,
                                         max_iter, points, weights);

      // renormalize if needed
      const double weights_sum = accumulate(begin(weights), end(weights), 0.0);
      if (weights_sum != 1.0)
        for (size_t i = 0; i < weights.size(); i++)
          weights[i] = weights[i]/weights_sum;

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
        estimated_unobs += counts_hist[1]*weights[i]/points[i];

      if (estimated_unobs > 0.0)
        estimated_unobs += distinct_obs;
      else {
        estimated_unobs = distinct_obs;
        n_points = 0;
      }

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "quadrature_estimated_unobs" << '\t' << "n_points" << endl;
      out << estimated_unobs << '\t' << n_points << endl;

    }
    // NOT QUICK MODE, BOOTSTRAP
    else {

      vector<double> quad_estimates;

      //setup rng
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

      for (size_t iter = 0; iter < max_iter && quad_estimates.size() < n_bootstraps; ++iter) {
        if (VERBOSE)
          cerr << "iter=" << "\t" << iter << endl;

        vector<double> sample_hist;
        resample_hist(rng, counts_hist_distinct_counts,
                      distinct_counts_hist, sample_hist);

        const double sampled_distinct = accumulate(begin(sample_hist), end(sample_hist), 0.0);

        // initialize moments, 0th moment is 1
        vector<double> bootstrap_moments(1, 1.0);
        // moments[r] = (r + 1)! n_{r+1} / n_1
        for (size_t i = 0; i < 2*max_num_points; i++) {
          bootstrap_moments.push_back(exp(factorial(i + 3) +
                                          log(sample_hist[i + 2]) -
                                          log(sample_hist[1])) );
        }

        size_t n_points = 0;
        n_points = ensure_pos_def_mom_seq(bootstrap_moments, tolerance, VERBOSE);
        n_points = min(n_points, max_num_points);
        if (VERBOSE)
          cerr << "n_points = " << n_points << endl;

        MomentSequence bootstrap_mom_seq(bootstrap_moments);

        vector<double> points, weights;
        bootstrap_mom_seq.Lower_quadrature_rules(VERBOSE, n_points, tolerance,
                                                 max_iter, points, weights);

        // renormalize if needed
        const double weights_sum = accumulate(begin(weights), end(weights), 0.0);
        if (weights_sum != 1.0)
          for (size_t i = 0; i < weights.size(); i++)
            weights[i] = weights[i]/weights_sum;

        double estimated_unobs = 0.0;

        for (size_t i = 0; i < weights.size(); i++)
          estimated_unobs += counts_hist[1]*weights[i]/points[i];

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
      median_and_ci(quad_estimates, c_level, median_estimate,
                    lower_ci, upper_ci);

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "median_estimated_unobs" << '\t'
          << "lower_ci" << '\t' << "upper_ci" << endl;
      out << median_estimate << '\t'
          << lower_ci << '\t' << upper_ci << endl;
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

int
main(const int argc, const char **argv) {

  static const string usage_message =
    "preseq:  a program for analyzing library complexity\n"
    "Version: " + preseq_version + "\n\n"
    "Usage:   preseq <command> [OPTIONS]\n\n"
    "<command>: c_curve    generate complexity curve for a library\n"
    "           lc_extrap  predict the yield for future experiments\n"
    "           gc_extrap  predict genome coverage low input\n"
    "                      sequencing experiments\n"
    "           bound_pop  lower bound on population size\n"
    "           pop_size   estimate number of unique species\n";

  if (argc < 2)
    cerr << usage_message << endl;

  else if (strcmp(argv[1], "lc_extrap") == 0) {

    return lc_extrap(0, argc, argv);

  }
  else if (strcmp(argv[1], "c_curve") == 0) {

    return c_curve(argc, argv);

  }
  else if (strcmp(argv[1], "gc_extrap") == 0) {

    return gc_extrap(argc, argv);

  }
  else if (strcmp(argv[1], "bound_pop") == 0) {

    return bound_pop(argc, argv);

  }
  else if (strcmp(argv[1], "pop_size") == 0) {

    return lc_extrap(1, argc, argv);

  }
  else {
    cerr << "unrecognized command: " << argv[1] << endl
         << usage_message << endl;
    return EXIT_SUCCESS;

  }
}
