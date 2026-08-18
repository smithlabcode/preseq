/* Copyright (C) 2013-2026 Andrew D. Smith and Timothy Daley
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "continued_fraction.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

#include "CLI11/CLI11.hpp"

#include <config.h>

#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <numeric>
#include <print>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

static constexpr auto about_msg = R"(
Predict properties of sequencing libraries.
)";

[[nodiscard]] static inline auto
rlstrip(const std::string &s) noexcept -> std::string {
  constexpr auto is_graph = [](const auto c) { return std::isgraph(c); };
  const auto start_itr = std::ranges::find_if(s, is_graph);
  auto stop_itr = std::end(s);
  while (stop_itr != std::cbegin(s) && !is_graph(*(stop_itr - 1)))
    --stop_itr;
  return std::string(start_itr, stop_itr);
}

template <typename T>
auto
get_counts_from_hist(const std::vector<T> &h) -> T {
  T c = 0.0;
  for (std::size_t i = 0; i < std::size(h); ++i)
    c += i * h[i];
  return c;
}

template <typename T>
auto
median_from_sorted_vector(const std::vector<T> sorted_data,
                          const std::size_t n) -> T {
  if (n == 0 || sorted_data.empty())
    return 0.0;

  const std::size_t lhs = (n - 1) / 2;
  const std::size_t rhs = n / 2;

  if (lhs == rhs)
    return sorted_data[lhs];

  return (sorted_data[lhs] + sorted_data[rhs]) / 2.0;
}

template <typename T>
auto
quantile_from_sorted_vector(const std::vector<T> sorted_data,
                            const std::size_t n, const double f) -> T {
  const double index = f * (n - 1);
  const std::size_t lhs = index;
  const double delta = index - lhs;

  if (n == 0 || sorted_data.empty())
    return 0.0;

  if (lhs + 1 == n)
    return sorted_data[lhs];

  return (1.0 - delta) * sorted_data[lhs] + delta * sorted_data[lhs + 1];
}

// Confidence interval stuff
static void
median_and_ci(std::vector<double> estimates,  // by val so we can sort them
              const double ci_level, double &median_estimate,
              double &lower_ci_estimate, double &upper_ci_estimate) {
  assert(std::size(estimates) > 0);
  std::ranges::sort(estimates);
  const double alpha = 1.0 - ci_level;
  const std::size_t N = std::size(estimates);
  median_estimate = median_from_sorted_vector(estimates, N);
  lower_ci_estimate = quantile_from_sorted_vector(estimates, N, alpha / 2);
  upper_ci_estimate =
    quantile_from_sorted_vector(estimates, N, 1.0 - alpha / 2);
}

static void
vector_median_and_ci(
  const std::vector<std::vector<double>> &bootstrap_estimates,
  const double ci_level, std::vector<double> &yield_estimates,
  std::vector<double> &lower_ci_lognorm,
  std::vector<double> &upper_ci_lognorm) {
  yield_estimates.clear();
  lower_ci_lognorm.clear();
  upper_ci_lognorm.clear();
  assert(!bootstrap_estimates.empty());
  const std::size_t n_est = std::size(bootstrap_estimates);
  std::vector<double> estimates_row(n_est, 0.0);
  for (std::size_t i = 0; i < std::size(bootstrap_estimates[0]); i++) {
    // estimates is in wrong order, work locally on const val
    for (std::size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];
    double median_estimate{};
    double lower_ci_estimate{};
    double upper_ci_estimate{};
    median_and_ci(estimates_row, ci_level, median_estimate, lower_ci_estimate,
                  upper_ci_estimate);
    std::ranges::sort(estimates_row);
    yield_estimates.push_back(median_estimate);
    lower_ci_lognorm.push_back(lower_ci_estimate);
    upper_ci_lognorm.push_back(upper_ci_estimate);
  }
}

template <typename uint_type>
void
multinomial(std::mt19937 &gen, const std::vector<double> &mult_probs,
            uint_type trials, std::vector<uint_type> &result) {
  using binom_dist = std::binomial_distribution<std::uint32_t>;
  result.clear();
  result.resize(std::size(mult_probs));
  auto remaining_prob =
    std::accumulate(std::cbegin(mult_probs), std::cend(mult_probs), 0.0);
  auto r = std::begin(result);
  auto p = std::begin(mult_probs);
  while (p != std::cend(mult_probs)) {  // iterate to sample for each category
    *r = binom_dist(trials, (*p) / remaining_prob)(gen);  // take the sample
    remaining_prob -= *p;  // update remaining probability mass
    trials -= *r;          // update remaining trials needed
    ++p;
    ++r;
  }
  if (trials > 0)
    throw std::runtime_error("multinomial sampling failed");
}

// Lanczos approximation for gamma function for x >= 0.5 - essentially an
// approximation for (x-1)!
static auto
factorial(double x) -> double {
  // constants
  static constexpr double LogRootTwoPi = 0.9189385332046727;
  static constexpr double Euler = std::numbers::e;

  // Approximation for factorial is actually x-1
  x -= 1.0;

  // clang-format off
  const auto lanczos = std::vector<double>{
    0.99999999999980993227684700473478,
    676.520368121885098567009190444019,
    -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
    12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7,
  };
  // clang-format on
  double Ag = lanczos[0];

  for (std::size_t k = 1; k < std::size(lanczos); k++)
    Ag += lanczos[k] / (x + k);

  const double term1 = (x + 0.5) * std::log((x + 7.5) / Euler);
  const double term2 = LogRootTwoPi + std::log(Ag);

  return term1 + (term2 - 7.0);
}

//  Extrap mode below here

// vals_hist[j] = n_{j} = # (counts = j)
// vals_hist_distinct_counts[k] = kth index j s.t. vals_hist[j] > 0
// stores kth index of vals_hist that is positive
// distinct_counts_hist[k] = vals_hist[vals_hist_distinct_counts[k]]
// stores the kth positive value of vals_hist
static void
resample_hist(std::mt19937 &gen,
              const std::vector<std::size_t> &vals_hist_distinct_counts,
              const std::vector<double> &distinct_counts_hist,
              std::vector<double> &out_hist) {
  const std::size_t hist_size = std::size(distinct_counts_hist);
  std::vector<std::uint32_t> sample_distinct_counts_hist(hist_size, 0);

  std::uint32_t distinct = std::accumulate(std::cbegin(distinct_counts_hist),
                                           end(distinct_counts_hist), 0.0);

  multinomial(gen, distinct_counts_hist, distinct, sample_distinct_counts_hist);

  out_hist.clear();
  out_hist.resize(vals_hist_distinct_counts.back() + 1, 0.0);
  for (std::size_t i = 0; i < hist_size; i++)
    out_hist[vals_hist_distinct_counts[i]] = sample_distinct_counts_hist[i];
}

// interpolate by explicit calculating the expectation for sampling without
// replacement; see K.L Heck 1975
//
// N total sample size
// S the total number of distincts
// n sub sample size
static auto
interpolate_distinct(const std::vector<double> &hist, const std::size_t N,
                     const std::size_t S, const std::size_t n) -> double {
  const double denom =
    factorial(N + 1) - factorial(n + 1) - factorial(N - n + 1);
  std::vector<double> numer(std::size(hist), 0);
  for (std::size_t i = 1; i < std::size(hist); i++) {
    // N - i -n + 1 should be greater than 0
    if (N < i + n) {
      numer[i] = 0;
    }
    else {
      const double x =
        factorial(N - i + 1) - factorial(n + 1) - factorial(N - i - n + 1);
      numer[i] = std::exp(x - denom) * hist[i];
    }
  }
  return S - std::accumulate(std::cbegin(numer), std::cend(numer), 0);
}

static void
extrapolate_curve(const ContinuedFraction &the_cf,
                  const double initial_distinct, const double vals_sum,
                  const double initial_sample_size, const double step_size,
                  const double max_sample_size,
                  std::vector<double> &estimates) {
  double curr_samp_sz = initial_sample_size;
  while (curr_samp_sz < max_sample_size) {
    const double fold = (curr_samp_sz - vals_sum) / vals_sum;
    assert(fold >= 0.0);
    estimates.push_back(initial_distinct + fold * the_cf(fold));
    curr_samp_sz += step_size;
  }
}

[[nodiscard]] static auto
extrap_bootstrap(
  const bool VERBOSE, const bool allow_defects, const std::uint64_t seed,
  const std::vector<double> &orig_hist, const std::size_t n_bootstraps,
  const std::size_t orig_max_terms, const int diagonal,
  const double bin_step_size, const double max_extrap,
  const std::size_t max_iter) -> std::vector<std::vector<double>> {
  // setup rng
  srand(time(nullptr) + getpid());
  std::mt19937 rng(seed);

  // const double vals_sum = get_counts_from_hist(orig_hist);
  const double initial_distinct =
    std::accumulate(std::cbegin(orig_hist), std::cend(orig_hist), 0.0);

  std::vector<std::size_t> orig_hist_distinct_counts;
  std::vector<double> distinct_orig_hist;
  for (std::size_t i = 0; i < std::size(orig_hist); i++)
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }

  std::vector<std::vector<double>> bootstrap_estimates;
  for (std::size_t iter = 0;
       (iter < max_iter && std::size(bootstrap_estimates) < n_bootstraps);
       ++iter) {
    if (VERBOSE && iter > 0 && iter % 72 == 0)
      std::println(std::cerr);  // bootstrap success progress only 72 char wide

    std::vector<double> yield_vector;
    std::vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);

    const double sample_vals_sum = get_counts_from_hist(hist);

    // resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    // compute complexity curve by random sampling w/out replacement
    const std::size_t distinct =
      std::accumulate(std::cbegin(hist), std::cend(hist), 0.0);
    std::size_t curr_sample_sz = bin_step_size;
    while (curr_sample_sz < sample_vals_sum) {
      yield_vector.push_back(
        interpolate_distinct(hist, sample_vals_sum, distinct, curr_sample_sz));
      curr_sample_sz += bin_step_size;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::size_t first_zero = 1;
    while (first_zero < std::size(hist) && hist[first_zero] > 0)
      ++first_zero;

    std::size_t max_terms = std::min(orig_max_terms, first_zero - 1);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 1);

    bool successful_bootstrap = false;
    // defect mode, simple extrapolation
    if (allow_defects) {
      const auto defect_cf = ContinuedFraction(hist, max_terms);
      defect_cf.extrapolate_curve(initial_distinct, sample_vals_sum,
                                  curr_sample_sz, bin_step_size, max_extrap,
                                  yield_vector);
      // no checking of curve in defect mode
      bootstrap_estimates.push_back(yield_vector);
      successful_bootstrap = true;
    }
    else {
      // refit curve for lower bound
      const ContinuedFractionApproximation lower_cfa(diagonal, max_terms);
      const auto lower_cf = lower_cfa.optimal_cf_distinct(hist);
      // extrapolate the curve start
      if (lower_cf.is_valid()) {
        lower_cf.extrapolate_curve(initial_distinct, sample_vals_sum,
                                   curr_sample_sz, bin_step_size, max_extrap,
                                   yield_vector);
        // sanity check
        if (check_yield_estimates_stability(yield_vector)) {
          bootstrap_estimates.push_back(yield_vector);
          successful_bootstrap = true;
        }
      }
    }
    if (VERBOSE)
      std::println(std::cerr, "{}", successful_bootstrap ? '.' : '_');
  }
  if (VERBOSE)
    std::println(std::cerr);
  if (std::size(bootstrap_estimates) < n_bootstraps)
    throw std::runtime_error("too many defects in the approximation, "
                             "consider running in defect mode");
  return bootstrap_estimates;
}

static auto
extrap_single_estimate(const bool VERBOSE, const bool allow_defects,
                       std::vector<double> &hist, std::size_t max_terms,
                       const int diagonal, const double step_size,
                       const double max_extrap,
                       std::vector<double> &yield_estimate) -> bool {
  yield_estimate.clear();

  const double vals_sum = get_counts_from_hist(hist);
  const double initial_distinct =
    std::accumulate(std::cbegin(hist), std::cend(hist), 0.0);

  // interpolate complexity curve by random sampling w/out replacement
  auto upper_limit = static_cast<std::size_t>(vals_sum);
  auto step = static_cast<std::size_t>(step_size);
  auto sample = static_cast<std::size_t>(step_size);
  for (; sample < upper_limit; sample += step)
    yield_estimate.push_back(
      interpolate_distinct(hist, upper_limit, initial_distinct, sample));

  // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
  std::size_t first_zero = 1;
  while (first_zero < std::size(hist) && hist[first_zero] > 0)
    ++first_zero;

  // Ensure we are not using a zero term
  max_terms = std::min(max_terms, first_zero - 1);

  // refit curve for lower bound (degree of approx is 1 less than
  // max_terms)
  max_terms = max_terms - (max_terms % 2 == 1);

  if (allow_defects) {
    std::vector<double> ps_coeffs;
    for (std::size_t j = 1; j <= max_terms; j++)
      ps_coeffs.push_back(hist[j] * std::pow(-1.0, j + 1));

    const ContinuedFraction defect_cf(ps_coeffs, diagonal, max_terms);

    extrapolate_curve(defect_cf, initial_distinct, vals_sum, sample, step_size,
                      max_extrap, yield_estimate);

    if (VERBOSE)
      std::println(std::cerr, "{}", defect_cf);
    // NO FAIL! defect mode doesn't care about failure
  }
  else {
    const ContinuedFractionApproximation lower_cfa(diagonal, max_terms);
    const ContinuedFraction lower_cf(lower_cfa.optimal_cf_distinct(hist));

    // extrapolate curve
    if (lower_cf.is_valid()) {
      extrapolate_curve(lower_cf, initial_distinct, vals_sum, sample, step_size,
                        max_extrap, yield_estimate);
    }
    else {
      // FAIL!
      // lower_cf unacceptable, need to bootstrap to obtain estimates
      return false;
    }

    if (VERBOSE)
      std::println(std::cerr, "{}", lower_cf);
  }
  // SUCCESS!!
  return true;
}

static auto
GoodToulmin2xExtrap(const std::vector<double> &counts_hist) -> double {
  double two_fold_extrap = 0.0;
  for (std::size_t i = 0; i < std::size(counts_hist); i++)
    two_fold_extrap += pow(-1.0, i + 1) * counts_hist[i];
  return two_fold_extrap;
}

static void
write_predicted_complexity_curve(
  const std::string &outfile, const double c_level, const double step_size,
  const std::vector<double> &yield_estimates,
  const std::vector<double> &yield_lower_ci_lognorm,
  const std::vector<double> &yield_upper_ci_lognorm) {
  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOWER_" << c_level << "CI\t"
      << "UPPER_" << c_level << "CI" << std::endl;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << std::endl;
  for (std::size_t i = 0; i < std::size(yield_estimates); ++i)
    out << (i + 1) * step_size << '\t' << yield_estimates[i] << '\t'
        << yield_lower_ci_lognorm[i] << '\t' << yield_upper_ci_lognorm[i]
        << std::endl;
}

// ADS: functions same, header different (above and this one)
static void
write_predicted_coverage_curve(
  const std::string &outfile, const double c_level, const double base_step_size,
  const std::size_t bin_size, const std::vector<double> &cvrg_estimates,
  const std::vector<double> &cvrg_lower_ci_lognorm,
  const std::vector<double> &cvrg_upper_ci_lognorm) {
  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "TOTAL_BASES\tEXPECTED_COVERED_BASES\t"
      << "LOWER_" << 100 * c_level << "%CI\t"
      << "UPPER_" << 100 * c_level << "%CI" << std::endl;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << std::endl;
  for (std::size_t i = 0; i < std::size(cvrg_estimates); ++i)
    out << (i + 1) * base_step_size << '\t' << cvrg_estimates[i] * bin_size
        << '\t' << cvrg_lower_ci_lognorm[i] * bin_size << '\t'
        << cvrg_upper_ci_lognorm[i] * bin_size << std::endl;
}

static auto
lc_extrap(int argc, char *argv[]) -> int {
  try {
    static const std::size_t min_required_counts = 4;
    static const std::string min_required_counts_error_message =
      "max count before zero is less than min required count (" +
      std::to_string(min_required_counts) + ") duplicates removed";

    std::string outfile;
    std::string input_file_name;

    std::size_t orig_max_terms = 100;
    double max_extrap = 1.0e10;
    double step_size = 1e6;
    std::size_t n_bootstraps = 100;
    int diagonal = 0;
    double c_level = 0.95;
    std::uint64_t seed = 408;

    /* FLAGS */
    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool SINGLE_ESTIMATE = false;
    bool allow_defects = false;
#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    std::size_t MAX_SEGMENT_LENGTH = 5000;
#endif
    const auto description = R"(
Extrapolate the complexity of a library. This is the approach described in
Daley & Smith (2013). The method applies rational function approximation via
continued fractions with the original goal of estimating the number of
distinct reads that a sequencing library would yield upon deeper sequencing.
This method has been used for many different purposes since then.
)";
    CLI::App app{about_msg};
    argv = app.ensure_utf8(argv);
    // app.usage(usage);
    if (argc >= 2)
      app.footer(rlstrip(description));

    // clang-format off
    app.set_help_flag("-h,--help", "print a detailed help message and exit");
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "output filename (directory must exist)")
      ->option_text("FILE")
      ->required();
    app.add_option("-e,--extrap", max_extrap, "maximum extrapolation");
    app.add_option("-s,--step", step_size, "extrapolation step size");
    app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--cval", c_level, "level for confidence intervals");
    app.add_option("-x,--terms", orig_max_terms, "maximum terms in estimator");
    app.add_option("-r,--seed", seed, "seed for random number generator");
#ifdef HAVE_HTSLIB
    app.add_option("-B,--bam", BAM_FORMAT_INPUT, "input is in BAM format");
    app.add_option("-l,--seg_len", MAX_SEGMENT_LENGTH,
                   "maximum segment length when merging paired end bam reads");
#endif
    app.add_flag("-P,--pe", PAIRED_END, "input is paired end read file");
    app.add_flag("-V,--vals", VALS_INPUT,
                 "input is a text file containing only the observed counts");
    app.add_flag("-H,--hist", HIST_INPUT,
                   "input is a text file containing the observed histogram");
    app.add_flag("-Q,--quick", SINGLE_ESTIMATE,
                 "quick mode (no bootstraps) for confidence intervals");
    app.add_flag("-D,--defects", allow_defects, "no testing for defects");
    app.add_flag("-v,--verbose", VERBOSE, "print more info");
    // clang-format on

    if (argc < 2) {
      // std::println("{}", app.help());
      std::cout << app.help() << std::endl;
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    std::vector<double> counts_hist;
    std::size_t n_reads{};

    /************ loading input ***************************************/
    if (HIST_INPUT) {
      if (VERBOSE)
        std::cerr << "HIST_INPUT" << std::endl;
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (VERBOSE)
        std::cerr << "VALS_INPUT" << std::endl;
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BAM_INPUT" << std::endl;
      const std::size_t MAX_READS_TO_HOLD = 5000000;
      std::size_t n_paired = 0;
      std::size_t n_mates = 0;
      n_reads =
        load_counts_BAM_pe(input_file_name, MAX_SEGMENT_LENGTH,
                           MAX_READS_TO_HOLD, n_paired, n_mates, counts_hist);
      if (VERBOSE) {
        std::cerr << "MERGED PAIRED END READS = " << n_paired << std::endl;
        std::cerr << "MATES PROCESSED = " << n_mates << std::endl;
      }
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        std::cerr << "BAM_INPUT" << std::endl;
      n_reads = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BED_INPUT" << std::endl;
      n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (VERBOSE)
        std::cerr << "BED_INPUT" << std::endl;
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }
    /************ done loading input **********************************/

    const std::size_t max_observed_count = std::size(counts_hist) - 1;
    const double distinct_reads =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::size_t first_zero = 1;
    while (first_zero < std::size(counts_hist) && counts_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = std::min(orig_max_terms, first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

    const std::size_t distinct_counts =
      std::count_if(std::cbegin(counts_hist), std::cend(counts_hist),
                    [](const double x) { return x > 0.0; });

    if (VERBOSE)
      std::cerr << "TOTAL READS     = " << n_reads << std::endl
                << "DISTINCT READS  = " << distinct_reads << std::endl
                << "DISTINCT COUNTS = " << distinct_counts << std::endl
                << "MAX COUNT       = " << max_observed_count << std::endl
                << "COUNTS OF 1     = " << counts_hist[1] << std::endl
                << "MAX TERMS       = " << orig_max_terms << std::endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      std::cerr << "OBSERVED COUNTS (" << std::size(counts_hist) << ")"
                << std::endl;
      for (std::size_t i = 0; i < std::size(counts_hist); i++)
        if (counts_hist[i] > 0)
          std::cerr << i << '\t' << static_cast<std::size_t>(counts_hist[i])
                    << std::endl;
      std::cerr << std::endl;
    }

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(counts_hist);
    if (two_fold_extrap < 0.0)
      throw std::runtime_error(
        "Saturation expected at double initial sample size."
        " Unable to extrapolate");

    // const std::size_t total_reads = get_counts_from_hist(counts_hist);

    // assert(total_reads == n_reads); // ADS: why commented out?

    // check that min required count is satisfied
    if (orig_max_terms < min_required_counts)
      throw std::runtime_error(min_required_counts_error_message);

    if (VERBOSE)
      std::cerr << "[ESTIMATING YIELD CURVE]" << std::endl;
    std::vector<double> yield_estimates;

    if (SINGLE_ESTIMATE) {
      const bool single_estimate_success = extrap_single_estimate(
        VERBOSE, allow_defects, counts_hist, orig_max_terms, diagonal,
        step_size, max_extrap, yield_estimates);
      // IF FAILURE, EXIT
      if (!single_estimate_success)
        throw std::runtime_error("single estimate failed, run "
                                 "full mode for estimates");

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS\tEXPECTED_DISTINCT" << std::endl;
      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << std::endl;
      for (std::size_t i = 0; i < std::size(yield_estimates); ++i)
        out << (i + 1) * step_size << '\t' << yield_estimates[i] << std::endl;
    }
    else {
      if (VERBOSE)
        std::cerr << "[BOOTSTRAPPING HISTOGRAM]" << std::endl;

      const std::size_t max_iter = 100 * n_bootstraps;

      const auto bootstrap_estimates = extrap_bootstrap(
        VERBOSE, allow_defects, seed, counts_hist, n_bootstraps, orig_max_terms,
        diagonal, step_size, max_extrap, max_iter);

      if (VERBOSE)
        std::cerr << "[COMPUTING CONFIDENCE INTERVALS]" << std::endl;
      // yield ci
      std::vector<double> yield_upper_ci_lognorm, yield_lower_ci_lognorm;

      vector_median_and_ci(bootstrap_estimates, c_level, yield_estimates,
                           yield_lower_ci_lognorm, yield_upper_ci_lognorm);

      /////////////////////////////////////////////////////////////////////
      if (VERBOSE)
        std::cerr << "[WRITING OUTPUT]" << std::endl;

      write_predicted_complexity_curve(outfile, c_level, step_size,
                                       yield_estimates, yield_lower_ci_lognorm,
                                       yield_upper_ci_lognorm);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///// GC_EXTRAP: predicting genomic coverage
/////

static auto
gc_extrap(int argc, char *argv[]) -> int {
  try {
    const std::size_t MIN_REQUIRED_COUNTS = 4;

    int diagonal = 0;
    std::size_t orig_max_terms = 100;
    std::size_t bin_size = 10;
    bool VERBOSE = false;
    std::string outfile;
    std::string input_file_name;
    double base_step_size = 1.0e8;
    std::size_t max_width = 10000;
    bool SINGLE_ESTIMATE = false;
    double max_extrap = 1.0e12;
    std::size_t n_bootstraps = 100;
    std::uint64_t seed = 408;
    bool allow_defects = false;

    bool NO_SEQUENCE = false;
    double c_level = 0.95;

    constexpr auto description = R"(
Extrapolate the size of the covered genome by mapped reads. This approach is
described in Daley & Smith (2014). The method is the same as for lc_extrap:
using rational function approximation to a power-series expansion for the
number of "unobserved" bases in the initial sample. The gc_extrap method is
adapted to deal with individual nucleotides rather than distinct reads.
)";

    CLI::App app{about_msg};
    argv = app.ensure_utf8(argv);
    // app.usage(usage);
    if (argc >= 2)
      app.footer(description);

    // clang-format off
    app.set_help_flag("-h,--help", "print a detailed help message and exit");
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "coverage yield output file")
      ->option_text("FILE")
      ->required();
    app.add_option("-w,--max_width", max_width,
                   "max fragment length, set equal to read length for single end reads");
    app.add_option("-b,--bin_size", bin_size, "bin size");
    app.add_option("-e,--extrap", max_extrap, "maximum extrapolation in base pairs");
    app.add_option("-s,--step", base_step_size, "step size in bases between extrapolations");
    app.add_option("-n,--bootstraps", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--cval", c_level, "level for confidence intervals");
    app.add_option("-x,--terms", orig_max_terms, "maximum number of terms");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_flag("-B,--bed", NO_SEQUENCE, "input is in bed format without sequence information");
    app.add_flag("-Q,--quick", SINGLE_ESTIMATE,
                 "quick mode: run gc_extrap without bootstrapping for confidence intervals");
    app.add_flag("-D,--defects", allow_defects,
                 "defects mode to extrapolate without testing for defects");
    app.add_flag("-v,--verbose", VERBOSE, "print more info");
    // clang-format on

    if (argc < 2) {
      // std::println("{}", app.help());
      std::cout << app.help() << std::endl;
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    std::vector<double> coverage_hist;
    std::size_t n_reads = 0;
    if (VERBOSE)
      std::cerr << "LOADING READS" << std::endl;

    if (NO_SEQUENCE) {
      if (VERBOSE)
        std::cerr << "BED FORMAT" << std::endl;
      n_reads = load_coverage_counts_GR(input_file_name, seed, bin_size,
                                        max_width, coverage_hist);
    }
    else {
      if (VERBOSE)
        std::cerr << "MAPPED READ FORMAT" << std::endl;
      n_reads = load_coverage_counts_MR(input_file_name, seed, bin_size,
                                        max_width, coverage_hist);
    }

    const double total_bins = get_counts_from_hist(coverage_hist);

    const double distinct_bins =
      std::accumulate(coverage_hist.begin(), coverage_hist.end(), 0.0);

    const double avg_bins_per_read = total_bins / n_reads;
    double bin_step_size = base_step_size / bin_size;

    const std::size_t max_observed_count = std::size(coverage_hist) - 1;

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::size_t first_zero = 1;
    while (first_zero < std::size(coverage_hist) &&
           coverage_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = std::min(orig_max_terms, first_zero - 1);

    if (VERBOSE)
      std::cerr << "TOTAL READS         = " << n_reads << std::endl
                << "BASE STEP SIZE      = " << base_step_size << std::endl
                << "BIN STEP SIZE       = " << bin_step_size << std::endl
                << "TOTAL BINS          = " << total_bins << std::endl
                << "BINS PER READ       = " << avg_bins_per_read << std::endl
                << "DISTINCT BINS       = " << distinct_bins << std::endl
                << "TOTAL BASES         = " << total_bins * bin_size
                << std::endl
                << "TOTAL COVERED BASES = " << distinct_bins * bin_size
                << std::endl
                << "MAX COVERAGE COUNT  = " << max_observed_count << std::endl
                << "COUNTS OF 1         = " << coverage_hist[1] << std::endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      std::cerr << "OBSERVED BIN COUNTS (" << std::size(coverage_hist) << ")"
                << std::endl;
      for (std::size_t i = 0; i < std::size(coverage_hist); i++)
        if (coverage_hist[i] > 0)
          std::cerr << i << '\t' << coverage_hist[i] << std::endl;
      std::cerr << std::endl;
    }

    // catch if all reads are distinct
    if (orig_max_terms < MIN_REQUIRED_COUNTS)
      throw std::runtime_error("max count before zero is les than min required "
                               "count (4), sample not sufficiently deep or "
                               "duplicates removed");

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(coverage_hist);
    if (two_fold_extrap < 0.0)
      throw std::runtime_error("Library expected to saturate in doubling of "
                               "experiment size, unable to extrapolate");

    if (VERBOSE)
      std::cerr << "[ESTIMATING COVERAGE CURVE]" << std::endl;

    std::vector<double> coverage_estimates;

    if (SINGLE_ESTIMATE) {
      bool SINGLE_ESTIMATE_SUCCESS = extrap_single_estimate(
        VERBOSE, allow_defects, coverage_hist, orig_max_terms, diagonal,
        bin_step_size, max_extrap / bin_size, coverage_estimates);
      // IF FAILURE, EXIT
      if (!SINGLE_ESTIMATE_SUCCESS)
        throw std::runtime_error("SINGLE ESTIMATE FAILED, NEED TO RUN IN "
                                 "FULL MODE FOR ESTIMATES");

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_BASES\tEXPECTED_DISTINCT" << std::endl;

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << std::endl;
      for (std::size_t i = 0; i < std::size(coverage_estimates); ++i)
        out << (i + 1) * base_step_size << '\t'
            << coverage_estimates[i] * bin_size << std::endl;
    }
    else {
      if (VERBOSE)
        std::cerr << "[BOOTSTRAPPING HISTOGRAM]" << std::endl;

      const std::size_t max_iter = 10 * n_bootstraps;

      const auto bootstrap_estimates =
        extrap_bootstrap(VERBOSE, allow_defects, seed, coverage_hist,
                         n_bootstraps, orig_max_terms, diagonal, bin_step_size,
                         max_extrap / bin_size, max_iter);

      if (VERBOSE)
        std::cerr << "[COMPUTING CONFIDENCE INTERVALS]" << std::endl;
      std::vector<double> coverage_upper_ci_lognorm, coverage_lower_ci_lognorm;
      vector_median_and_ci(bootstrap_estimates, c_level, coverage_estimates,
                           coverage_lower_ci_lognorm,
                           coverage_upper_ci_lognorm);

      if (VERBOSE)
        std::cerr << "[WRITING OUTPUT]" << std::endl;

      write_predicted_coverage_curve(
        outfile, c_level, base_step_size, bin_size, coverage_estimates,
        coverage_lower_ci_lognorm, coverage_upper_ci_lognorm);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

static auto
c_curve(int argc, char *argv[]) -> int {
  try {
    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    std::uint64_t seed = 408;
    std::string outfile;
    std::string input_file_name;
    std::size_t upper_limit = 0;
    double step_size = 1e6;
#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    std::size_t MAX_SEGMENT_LENGTH = 5000;
#endif
    const auto description = R"(
Generate the complexity curve for data. This does not extrapolate, but instead
resamples from the given data.
)";
    CLI::App app{about_msg};
    argv = app.ensure_utf8(argv);
    // app.usage(usage);
    if (argc >= 2)
      app.footer(description);

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
    app.add_option("-v,--verbose", VERBOSE, "print more info");
    // clang-format on

    if (argc < 2) {
      // std::println("{}", app.help());
      std::cout << app.help() << std::endl;
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    // Setup the random number generator
    srand(time(nullptr) + getpid());  // give the random fxn a new seed
    std::mt19937 rng(seed);

    std::vector<double> counts_hist;
    std::size_t n_reads = 0;

    // LOAD VALUES
    if (HIST_INPUT) {
      if (VERBOSE)
        std::cerr << "INPUT_HIST" << std::endl;
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (VERBOSE)
        std::cerr << "VALS_INPUT" << std::endl;
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BAM_INPUT" << std::endl;
      const std::size_t MAX_READS_TO_HOLD = 5000000;
      std::size_t n_paired = 0;
      std::size_t n_mates = 0;
      n_reads =
        load_counts_BAM_pe(input_file_name, MAX_SEGMENT_LENGTH,
                           MAX_READS_TO_HOLD, n_paired, n_mates, counts_hist);
      if (VERBOSE)
        std::cerr << "MERGED PAIRED END READS = " << n_paired << std::endl
                  << "MATES PROCESSED = " << n_mates << std::endl;
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        std::cerr << "BAM_INPUT" << std::endl;
      n_reads = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BED_INPUT" << std::endl;
      n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (VERBOSE)
        std::cerr << "BED_INPUT" << std::endl;
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }

    const std::size_t max_observed_count = std::size(counts_hist) - 1;
    const double distinct_reads =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    const std::size_t total_reads = get_counts_from_hist(counts_hist);

    const std::size_t distinct_counts =
      std::count_if(std::cbegin(counts_hist), std::cend(counts_hist),
                    [](const double x) { return x > 0.0; });

    if (VERBOSE)
      std::cerr << "TOTAL READS     = " << n_reads << std::endl
                << "COUNTS_SUM      = " << total_reads << std::endl
                << "DISTINCT READS  = " << distinct_reads << std::endl
                << "DISTINCT COUNTS = " << distinct_counts << std::endl
                << "MAX COUNT       = " << max_observed_count << std::endl
                << "COUNTS OF 1     = " << counts_hist[1] << std::endl;

    if (VERBOSE) {
      // output the original histogram
      std::cerr << "OBSERVED COUNTS (" << std::size(counts_hist) << ")"
                << std::endl;
      for (std::size_t i = 0; i < std::size(counts_hist); i++)
        if (counts_hist[i] > 0)
          std::cerr << i << '\t' << static_cast<std::size_t>(counts_hist[i])
                    << std::endl;
      std::cerr << std::endl;
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
    out << "total_reads" << "\t" << "distinct_reads" << std::endl;
    out << 0 << '\t' << 0 << std::endl;
    for (std::size_t i = step_size; i <= upper_limit; i += step_size) {
      if (VERBOSE)
        std::cerr << "sample size: " << i << std::endl;
      out << i << "\t"
          << interpolate_distinct(counts_hist, total_reads, distinct_reads, i)
          << std::endl;
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// BOUND_UNOBS: bounding n_0

static auto
bound_pop(int argc, char *argv[]) -> int {
  try {
    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    bool QUICK_MODE = false;

    std::string outfile;
    std::string input_file_name;

#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    std::size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    std::size_t max_num_points = 10;
    double tolerance = 1e-20;
    std::size_t n_bootstraps = 500;
    double c_level = 0.95;
    std::size_t max_iter = 100;
    std::uint64_t seed = 408;

    const auto description = R"(
Estimate the size of the underlying population based on counts
of observed species in an initial sample.
)";
    CLI::App app{about_msg};
    argv = app.ensure_utf8(argv);
    // app.usage(usage);
    if (argc >= 2)
      app.footer(description);

    // clang-format off
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "species richness output file (default: stdout)");
    app.add_option("-p,--max_num_points", max_num_points, "maximum number of points in quadrature estimates");
    app.add_option("-t,--tolerance", tolerance, "numerical tolerance");
    app.add_option("-n,--bootstraps", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--clevel", c_level, "level for confidence intervals");
    app.add_option("-P,--pe", PAIRED_END, "input is paired end read file");
    app.add_option("-H,--hist", HIST_INPUT, "input is a text file containing the observed histogram");
    app.add_option("-V,--vals", VALS_INPUT, "input is a text file containing only the observed duplicate counts");
#ifdef HAVE_HTSLIB
    app.add_option("-B,--bam", BAM_FORMAT_INPUT, "input is in BAM format");
    app.add_option("-l,--seg_len", MAX_SEGMENT_LENGTH, "maximum segment length when merging paired end bam reads");
#endif
    app.add_option("-Q,--quick", QUICK_MODE, "quick mode, estimate without bootstrapping");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_option("-v,--verbose", VERBOSE, "print more info");
    // clang-format on

    if (argc < 2) {
      // std::println("{}", app.help());
      std::cout << app.help() << std::endl;
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    std::vector<double> counts_hist;
    std::size_t n_obs = 0;

    // LOAD VALUES
    if (HIST_INPUT) {
      if (VERBOSE)
        std::cerr << "HIST_INPUT" << std::endl;
      n_obs = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (VERBOSE)
        std::cerr << "VALS_INPUT" << std::endl;
      n_obs = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BAM_INPUT" << std::endl;
      const std::size_t MAX_READS_TO_HOLD = 5000000;
      std::size_t n_paired = 0;
      std::size_t n_mates = 0;
      n_obs =
        load_counts_BAM_pe(input_file_name, MAX_SEGMENT_LENGTH,
                           MAX_READS_TO_HOLD, n_paired, n_mates, counts_hist);
      if (VERBOSE) {
        std::cerr << "MERGED PAIRED END READS = " << n_paired << std::endl;
        std::cerr << "MATES PROCESSED = " << n_mates << std::endl;
      }
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        std::cerr << "BAM_INPUT" << std::endl;
      n_obs = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BED_INPUT" << std::endl;
      n_obs = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (VERBOSE)
        std::cerr << "BED_INPUT" << std::endl;
      n_obs = load_counts_BED_se(input_file_name, counts_hist);
    }

    const double distinct_obs =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    std::vector<double> measure_moments;
    // mu_r = (r + 1)! n_{r+1} / n_1
    std::size_t idx = 1;
    while (counts_hist[idx] > 0 && idx <= std::size(counts_hist)) {
      // idx + 1 because function calculates (x-1)!
      measure_moments.push_back(std::exp(
        factorial(idx + 1) + log(counts_hist[idx]) - log(counts_hist[1])));
      if (!std::isfinite(measure_moments.back())) {
        measure_moments.pop_back();
        break;
      }
      idx++;
    }

    if (VERBOSE) {
      std::cerr << "TOTAL OBSERVATIONS     = " << n_obs << std::endl
                << "DISTINCT OBSERVATIONS  = " << distinct_obs << std::endl
                << "MAX COUNT              = " << std::size(counts_hist) - 1
                << std::endl;

      // OUTPUT THE ORIGINAL HISTOGRAM
      std::cerr << "OBSERVED COUNTS (" << std::size(counts_hist) << ")"
                << std::endl;
      for (std::size_t i = 0; i < std::size(counts_hist); i++)
        if (counts_hist[i] > 0)
          std::cerr << i << '\t' << std::setprecision(16) << counts_hist[i]
                    << std::endl;

      std::cerr << "OBSERVED MOMENTS" << std::endl;
      for (const double measure_moment : measure_moments)
        std::cerr << std::setprecision(16) << measure_moment << std::endl;
    }

    if (QUICK_MODE) {
      if (std::size(measure_moments) < 2 * max_num_points)
        max_num_points =
          static_cast<std::size_t>(floor(std::size(measure_moments) / 2));
      else
        measure_moments.resize(2 * max_num_points);
      std::size_t n_points = 0;
      n_points = ensure_pos_def_mom_seq(measure_moments, tolerance, VERBOSE);
      if (VERBOSE)
        std::cerr << "n_points = " << n_points << std::endl;

      MomentSequence obs_mom_seq(measure_moments);

      if (VERBOSE) {
        for (std::size_t k = 0; k < std::size(obs_mom_seq.alpha); k++)
          std::cerr << "alpha_" << k << '\t';
        std::cerr << std::endl;
        for (const double k : obs_mom_seq.alpha)
          std::cerr << k << '\t';
        std::cerr << std::endl;

        for (std::size_t k = 0; k < std::size(obs_mom_seq.beta); k++)
          std::cerr << "beta_" << k << '\t';
        std::cerr << std::endl;
        for (const double k : obs_mom_seq.beta)
          std::cerr << k << '\t';
        std::cerr << std::endl;
      }

      std::vector<double> points, weights;
      obs_mom_seq.Lower_quadrature_rules(n_points, tolerance, max_iter, points,
                                         weights);

      // renormalize if needed
      const double weights_sum =
        std::accumulate(std::cbegin(weights), std::cend(weights), 0.0);
      if (weights_sum != 1.0)
        for (double &weight : weights)
          weight = weight / weights_sum;

      if (VERBOSE) {
        std::cerr << "points = " << std::endl;
        for (const double point : points)
          std::cerr << point << '\t';
        std::cerr << std::endl;

        std::cerr << "weights = " << std::endl;
        for (const double weight : weights)
          std::cerr << weight << '\t';
        std::cerr << std::endl;
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
        of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "quadrature_estimated_unobs" << '\t' << "n_points" << std::endl;
      out << estimated_unobs << '\t' << n_points << std::endl;
    }
    // NOT QUICK MODE, BOOTSTRAP
    else {
      std::vector<double> quad_estimates;

      // setup rng
      srand(time(nullptr) + getpid());
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
        if (VERBOSE)
          std::cerr << "iter=" << "\t" << iter << std::endl;

        std::vector<double> sample_hist;
        resample_hist(rng, counts_hist_distinct_counts, distinct_counts_hist,
                      sample_hist);

        const double sampled_distinct = std::accumulate(
          std::cbegin(sample_hist), std::cend(sample_hist), 0.0);

        // initialize moments, 0th moment is 1
        std::vector<double> bootstrap_moments(1, 1.0);
        // moments[r] = (r + 1)! n_{r+1} / n_1
        for (std::size_t i = 0; i < 2 * max_num_points; i++) {
          bootstrap_moments.push_back(std::exp(
            factorial(i + 3) + log(sample_hist[i + 2]) - log(sample_hist[1])));
        }

        std::size_t n_points = 0;
        n_points =
          ensure_pos_def_mom_seq(bootstrap_moments, tolerance, VERBOSE);
        n_points = std::min(n_points, max_num_points);
        if (VERBOSE)
          std::cerr << "n_points = " << n_points << std::endl;

        MomentSequence bootstrap_mom_seq(bootstrap_moments);

        std::vector<double> points, weights;
        bootstrap_mom_seq.Lower_quadrature_rules(n_points, tolerance, max_iter,
                                                 points, weights);

        // renormalize if needed
        const double weights_sum =
          std::accumulate(std::cbegin(weights), std::cend(weights), 0.0);
        if (weights_sum != 1.0)
          for (double &weight : weights)
            weight = weight / weights_sum;

        double estimated_unobs = 0.0;

        for (std::size_t i = 0; i < std::size(weights); i++)
          estimated_unobs += counts_hist[1] * weights[i] / points[i];

        if (estimated_unobs > 0.0)
          estimated_unobs += sampled_distinct;
        else {
          estimated_unobs = sampled_distinct;
          n_points = 0;
        }

        if (VERBOSE) {
          std::cerr << "bootstrapped_moments=" << std::endl;
          for (const double bootstrap_moment : bootstrap_moments)
            std::cerr << bootstrap_moment << std::endl;
        }
        if (VERBOSE) {
          for (std::size_t k = 0; k < std::size(bootstrap_mom_seq.alpha); k++)
            std::cerr << "alpha_" << k << '\t';
          std::cerr << std::endl;
          for (const double k : bootstrap_mom_seq.alpha)
            std::cerr << k << '\t';
          std::cerr << std::endl;

          for (std::size_t k = 0; k < std::size(bootstrap_mom_seq.beta); k++)
            std::cerr << "beta_" << k << '\t';
          std::cerr << std::endl;
          for (const double k : bootstrap_mom_seq.beta)
            std::cerr << k << '\t';
          std::cerr << std::endl;
        }
        if (VERBOSE) {
          std::cerr << "points=" << "\t";
          for (const double point : points)
            std::cerr << point << "\t";
          std::cerr << std::endl;
          std::cerr << "weights=" << "\t";
          for (const double weight : weights)
            std::cerr << weight << "\t";
          std::cerr << std::endl;
          std::cerr << "estimated_unobs=" << "\t" << estimated_unobs
                    << std::endl;
        }

        quad_estimates.push_back(estimated_unobs);
      }

      double median_estimate, lower_ci, upper_ci;
      median_and_ci(quad_estimates, c_level, median_estimate, lower_ci,
                    upper_ci);

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "median_estimated_unobs" << '\t' << "lower_ci" << '\t'
          << "upper_ci" << std::endl;
      out << median_estimate << '\t' << lower_ci << '\t' << upper_ci
          << std::endl;
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

static auto
pop_size(int argc, char *argv[]) -> int {
  try {
    static const std::size_t min_required_counts = 4;
    static const std::string min_required_counts_error_message =
      "max count before zero is less than min required count (" +
      std::to_string(min_required_counts) + ") duplicates removed";

    std::string outfile;
    std::string input_file_name;

    std::size_t orig_max_terms = 100;
    double max_extrap = 0.0;
    double step_size = 0.0;
    std::size_t n_desired_steps = 50;
    std::size_t n_bootstraps = 100;
    int diagonal = 0;
    double c_level = 0.95;
    std::uint64_t seed = 408;

    /* FLAGS */
    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool SINGLE_ESTIMATE = false;
    bool allow_defects = false;

#ifdef HAVE_HTSLIB
    bool BAM_FORMAT_INPUT = false;
    std::size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    const auto description = R"(
Estimate the total population size using the approach described in
Daley & Smith (2013), extrapolating to very long range. Default
parameters assume that the initial sample represents at least
1e-9 of the population, which is sufficient for every example
application we have seen.
)";
    CLI::App app{about_msg};
    argv = app.ensure_utf8(argv);
    // app.usage(usage);
    if (argc >= 2)
      app.footer(description);

    // clang-format off
    app.add_option("-i,--input", input_file_name, "input file")
      ->option_text("FILE")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "yield output file default: stdout");
    app.add_option("-e,--extrap", max_extrap, "maximum extrapolation");
    app.add_option("-s,--steps", n_desired_steps, "number of steps");
    app.add_option("-n,--boots", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--cval", c_level, "level for confidence intervals");
    app.add_option("-x,--terms", orig_max_terms, "maximum terms in estimator");
#ifdef HAVE_HTSLIB
    app.add_option("-B,--bam", BAM_FORMAT_INPUT, "input is in BAM format");
    app.add_option("-l,--seg_len", MAX_SEGMENT_LENGTH, "maximum segment length when merging paired end bam reads");
#endif
    app.add_option("-P,--pe", PAIRED_END, "input is paired end read file");
    app.add_option("-V,--vals", VALS_INPUT, "input is a text file containing only the observed counts");
    app.add_option("-H,--hist", HIST_INPUT, "input is a text file containing the observed histogram");
    app.add_option("-Q,--quick", SINGLE_ESTIMATE,
                   "quick mode (no bootstraps) for confidence intervals");
    app.add_option("-D,--defects", allow_defects, "no testing for defects");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_option("-v,--verbose", VERBOSE, "print more info");
    // clang-format on

    if (argc < 2) {
      // std::println("{}", app.help());
      std::cout << app.help() << std::endl;
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    std::vector<double> counts_hist;
    std::size_t n_reads = 0;

    /************ loading input ***************************************/
    if (HIST_INPUT) {
      if (VERBOSE)
        std::cerr << "HIST_INPUT" << std::endl;
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (VERBOSE)
        std::cerr << "VALS_INPUT" << std::endl;
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_HTSLIB
    else if (BAM_FORMAT_INPUT && PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BAM_INPUT" << std::endl;
      const std::size_t MAX_READS_TO_HOLD = 5000000;
      std::size_t n_paired = 0;
      std::size_t n_mates = 0;
      n_reads =
        load_counts_BAM_pe(input_file_name, MAX_SEGMENT_LENGTH,
                           MAX_READS_TO_HOLD, n_paired, n_mates, counts_hist);
      if (VERBOSE) {
        std::cerr << "MERGED PAIRED END READS = " << n_paired << std::endl;
        std::cerr << "MATES PROCESSED = " << n_mates << std::endl;
      }
    }
    else if (BAM_FORMAT_INPUT) {
      if (VERBOSE)
        std::cerr << "BAM_INPUT" << std::endl;
      n_reads = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if (PAIRED_END) {
      if (VERBOSE)
        std::cerr << "PAIRED_END_BED_INPUT" << std::endl;
      n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else {  // default is single end bed file
      if (VERBOSE)
        std::cerr << "BED_INPUT" << std::endl;
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }
    /************ done loading input **********************************/

    const std::size_t max_observed_count = std::size(counts_hist) - 1;
    const double distinct_reads =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    std::size_t first_zero = 1;
    while (first_zero < std::size(counts_hist) && counts_hist[first_zero] > 0)
      ++first_zero;

    orig_max_terms = std::min(orig_max_terms, first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);

    if (max_extrap < 1.0)
      max_extrap = 1000000000 * distinct_reads;
    if (step_size < 1.0)
      step_size = (max_extrap - distinct_reads) / n_desired_steps;

    const std::size_t distinct_counts =
      std::count_if(std::cbegin(counts_hist), std::cend(counts_hist),
                    [](const double x) { return x > 0.0; });

    if (VERBOSE)
      std::cerr << "TOTAL READS     = " << n_reads << std::endl
                << "DISTINCT READS  = " << distinct_reads << std::endl
                << "DISTINCT COUNTS = " << distinct_counts << std::endl
                << "MAX COUNT       = " << max_observed_count << std::endl
                << "COUNTS OF 1     = " << counts_hist[1] << std::endl
                << "MAX TERMS       = " << orig_max_terms << std::endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      std::cerr << "OBSERVED COUNTS (" << std::size(counts_hist) << ")"
                << std::endl;
      for (std::size_t i = 0; i < std::size(counts_hist); i++)
        if (counts_hist[i] > 0)
          std::cerr << i << '\t' << static_cast<std::size_t>(counts_hist[i])
                    << std::endl;
      std::cerr << std::endl;
    }

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(counts_hist);
    if (two_fold_extrap < 0.0)
      throw std::runtime_error(
        "Saturation expected at double initial sample size."
        " Unable to extrapolate");

    // const std::size_t total_reads = get_counts_from_hist(counts_hist);

    // assert(total_reads == n_reads); // ADS: why commented out?

    // check that min required count is satisfied
    if (orig_max_terms < min_required_counts)
      throw std::runtime_error(min_required_counts_error_message);

    if (VERBOSE)
      std::cerr << "[ESTIMATING YIELD CURVE]" << std::endl;

    std::vector<double> yield_estimates;

    if (SINGLE_ESTIMATE) {
      const bool single_estimate_success = extrap_single_estimate(
        VERBOSE, allow_defects, counts_hist, orig_max_terms, diagonal,
        step_size, max_extrap, yield_estimates);
      // IF FAILURE, EXIT
      if (!single_estimate_success)
        throw std::runtime_error("single estimate failed, run "
                                 "full mode for estimates");

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS\tEXPECTED_DISTINCT" << std::endl;
      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << std::endl;
      for (std::size_t i = 0; i < std::size(yield_estimates); ++i)
        out << (i + 1) * step_size << '\t' << yield_estimates[i] << std::endl;
    }
    else {
      if (VERBOSE)
        std::cerr << "[BOOTSTRAPPING HISTOGRAM]" << std::endl;

      const std::size_t max_iter = 100 * n_bootstraps;

      const auto bootstrap_estimates = extrap_bootstrap(
        VERBOSE, allow_defects, seed, counts_hist, n_bootstraps, orig_max_terms,
        diagonal, step_size, max_extrap, max_iter);

      if (VERBOSE)
        std::cerr << "[COMPUTING CONFIDENCE INTERVALS]" << std::endl;
      // yield ci
      std::vector<double> yield_upper_ci_lognorm, yield_lower_ci_lognorm;

      vector_median_and_ci(bootstrap_estimates, c_level, yield_estimates,
                           yield_lower_ci_lognorm, yield_upper_ci_lognorm);

      if (VERBOSE)
        std::cerr << "[WRITING OUTPUT]" << std::endl;

      std::ofstream of;
      if (!outfile.empty())
        of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      const std::size_t n_ests = std::size(yield_estimates) - 1;
      if (n_ests < 2)
        throw std::runtime_error(
          "problem with number of estimates in pop_size");

      const bool converged =
        (yield_estimates[n_ests] - yield_estimates[n_ests - 1] < 1.0);

      out << "pop_size_estimate" << '\t' << "lower_ci" << '\t' << "upper_ci"
          << std::endl;
      out << yield_estimates.back() << '\t' << yield_lower_ci_lognorm.back()
          << '\t' << yield_upper_ci_lognorm.back();
      if (!converged)
        out << "\tnot_converged";
      out << std::endl;
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int
main(int argc, char *argv[]) {
  static const auto usage_message = R"(
preseq: a program for analyzing library complexity

Usage: preseq <command> [OPTIONS]

<command>:  c_curve    generate complexity curve for a library
            lc_extrap  predict the yield for future experiments
            gc_extrap  predict genome coverage low input sequencing experiments
            bound_pop  lower bound on population size
            pop_size   estimate number of unique species

Version: )" + std::string(VERSION);

  if (argc < 2) {
    std::cerr << rlstrip(usage_message) << std::endl;
    return EXIT_SUCCESS;
  }

  if (std::strcmp(argv[1], "lc_extrap") == 0)
    return lc_extrap(argc - 1, argv + 1);

  if (std::strcmp(argv[1], "c_curve") == 0)
    return c_curve(argc - 1, argv + 1);

  if (std::strcmp(argv[1], "gc_extrap") == 0)
    return gc_extrap(argc - 1, argv + 1);

  if (std::strcmp(argv[1], "bound_pop") == 0)
    return bound_pop(argc - 1, argv + 1);

  if (std::strcmp(argv[1], "pop_size") == 0)
    return pop_size(argc - 1, argv + 1);

  std::cerr << "unrecognized command: " << argv[1] << std::endl
            << usage_message << std::endl;
  return EXIT_SUCCESS;
}
