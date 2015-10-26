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

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>
#include <queue>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <tr1/unordered_map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>

#define PRESEQ_VERSION "2.0.0"

// AS: might not be good to depend on mapped read here
// TD: if we're including gc_extrap, we need the dependence

#include "continued_fraction.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

using std::string;
using std::min;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::ifstream;
using std::isfinite;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;



/////////////////////////////////////////////////////////
// Confidence interval stuff

static inline double
alpha_log_confint_multiplier(const double estimate,
                             const double variance, const double alpha) {
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  return exp(inv_norm_alpha*
             sqrt(log(1.0 + variance/pow(estimate, 2))));
}


static void
median_and_ci(const vector<double> &estimates,
              const double ci_level,
              double &median_estimate,
              double &lower_ci_estimate,
              double &upper_ci_estimate){
  assert(!estimates.empty());
  const double alpha = 1.0 - ci_level;
  const size_t n_est = estimates.size();
  vector<double> sorted_estimates(estimates);
  sort(sorted_estimates.begin(), sorted_estimates.end());
  median_estimate =
    gsl_stats_median_from_sorted_data(&sorted_estimates[0], 
                                      1, n_est);
  const double variance = 
    gsl_stats_variance(&sorted_estimates[0], 1, n_est);
  const double confint_mltr =
    alpha_log_confint_multiplier(median_estimate, variance, alpha);

  lower_ci_estimate = median_estimate/confint_mltr;
  upper_ci_estimate = median_estimate*confint_mltr;

}

static void
vector_median_and_ci(const vector<vector<double> > &bootstrap_estimates,
                     const double ci_level, 
                     vector<double> &yield_estimates,
                     vector<double> &lower_ci_lognormal,
                     vector<double> &upper_ci_lognormal) {

  yield_estimates.clear();
  lower_ci_lognormal.clear();
  upper_ci_lognormal.clear();
  assert(!bootstrap_estimates.empty());

  const size_t n_est = bootstrap_estimates.size();
  vector<double> estimates_row(bootstrap_estimates.size(), 0.0);
  for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {

    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];

    double median_estimate, lower_ci_estimate, upper_ci_estimate;
    median_and_ci(estimates_row, ci_level, median_estimate,
                  lower_ci_estimate, upper_ci_estimate);
    sort(estimates_row.begin(), estimates_row.end());

    yield_estimates.push_back(median_estimate);
    lower_ci_lognormal.push_back(lower_ci_estimate);
    upper_ci_lognormal.push_back(upper_ci_estimate);
  }
}

void
log_mean(const bool VERBOSE,
	 const vector<double> &estimates,
	 const double c_level,
	 double &log_mean, 
	 double &log_lower_ci,
	 double &log_upper_ci){
  vector<double> log_estimates(estimates);
  for(size_t i = 0; i < log_estimates.size(); i++)
    log_estimates[i] = log(log_estimates[i]);

  log_mean = exp(gsl_stats_mean(&log_estimates[0], 1,
				log_estimates.size()) );

  double log_std_dev = std::sqrt(gsl_stats_variance(&log_estimates[0], 1, 
						    log_estimates.size()) );

  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv((1.0 - c_level)/2.0);
  log_lower_ci = exp(log(log_mean) - inv_norm_alpha*log_std_dev);
  log_upper_ci = exp(log(log_mean) + inv_norm_alpha*log_std_dev);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////  EXTRAP MODE BELOW HERE
/////


// vals_hist[j] = n_{j} = # (counts = j)
// vals_hist_distinct_counts[k] = kth index j s.t. vals_hist[j] > 0
// stores kth index of vals_hist that is positive
// distinct_counts_hist[k] = vals_hist[vals_hist_distinct_counts[k]]
// stores the kth positive value of vals_hist
void
resample_hist(const gsl_rng *rng, const vector<size_t> &vals_hist_distinct_counts,
              const vector<double> &distinct_counts_hist,
              vector<double> &out_hist) {

  vector<unsigned int> sample_distinct_counts_hist(distinct_counts_hist.size(), 0);

  const unsigned int distinct =
    static_cast<unsigned int>(accumulate(distinct_counts_hist.begin(),
                                         distinct_counts_hist.end(), 0.0));

  gsl_ran_multinomial(rng, distinct_counts_hist.size(), distinct,
                      &distinct_counts_hist.front(),
                      &sample_distinct_counts_hist.front());

  out_hist.clear();
  out_hist.resize(vals_hist_distinct_counts.back() + 1, 0.0);
  for(size_t i = 0; i < sample_distinct_counts_hist.size(); i++)
    out_hist[vals_hist_distinct_counts[i]] =
      static_cast<double>(sample_distinct_counts_hist[i]);
}

// interpolate by explicit calculating the expectation 
// for sampling without replacement; 
// see K.L Heck 1975
// N total sample size; S the total number of distincts
// n sub sample size
static double
interpolate_distinct(vector<double> &hist, size_t N,
                      size_t S, const size_t n) {
  double denom = gsl_sf_lngamma(N + 1) - gsl_sf_lngamma(n + 1) - gsl_sf_lngamma(N - n + 1);
  vector<double> numer(hist.size(), 0); 
  for (size_t i = 1; i < hist.size(); i++) {
	// N - i -n + 1 should be greater than 0
	if (N < i + n) {
	  numer[i] = 0;
	} else {
	  numer[i] = gsl_sf_lngamma(N - i + 1) - gsl_sf_lngamma(n + 1) - gsl_sf_lngamma(N - i - n + 1);
	  numer[i] = exp(numer[i] - denom) * hist[i];
	}
  }
  return S - accumulate(numer.begin(), numer.end(), 0);
}


// check if estimates are finite, increasing, and concave
static bool
check_yield_estimates(const vector<double> &estimates) {

  if (estimates.empty())
    return false;

  // make sure that the estimate is increasing in the time_step and is
  // below the initial distinct per step_size
  if (!isfinite(accumulate(estimates.begin(), estimates.end(), 0.0)))
    return false;

  for (size_t i = 1; i < estimates.size(); ++i)
    if ((estimates[i] < estimates[i - 1]) ||
        (i >= 2 && (estimates[i] - estimates[i - 1] >
                    estimates[i - 1] - estimates[i - 2])) ||
        (estimates[i] < 0.0))
      return false;

  return true;
}


void
extrap_bootstrap(const bool VERBOSE, const vector<double> &orig_hist,
                 const size_t bootstraps, const size_t orig_max_terms,
                 const int diagonal, const double bin_step_size,
                 const double max_extrapolation, const size_t max_iter,
                 vector< vector<double> > &bootstrap_estimates) {
  // clear returning vectors
  bootstrap_estimates.clear();

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand());

  double vals_sum = 0.0;
  for(size_t i = 0; i < orig_hist.size(); i++)
    vals_sum += orig_hist[i]*i;

  const double initial_distinct 
    = accumulate(orig_hist.begin(), orig_hist.end(), 0.0);


  vector<size_t> orig_hist_distinct_counts;
  vector<double> distinct_orig_hist;
  for (size_t i = 0; i < orig_hist.size(); i++){
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }
  }
  
  for (size_t iter = 0;
       (iter < max_iter && bootstrap_estimates.size() < bootstraps);
       ++iter) {

    vector<double> yield_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);

    double sample_vals_sum = 0.0;
    for(size_t i = 0; i < hist.size(); i++)
      sample_vals_sum += i*hist[i];

    //resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    // compute complexity curve by random sampling w/out replacement
    const size_t upper_limit = static_cast<size_t>(sample_vals_sum);
	const size_t distinct = static_cast<size_t>(accumulate(hist.begin(), hist.end(), 0.0));
    const size_t step = static_cast<size_t>(bin_step_size);
    size_t sample = step;
    while(sample < upper_limit){
      yield_vector.push_back(interpolate_distinct(hist, upper_limit, distinct, sample));
      sample += step;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < hist.size() &&
           hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
    
    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 1);
    
    //refit curve for lower bound
    const ContinuedFractionApproximation
      lower_cfa(diagonal, max_terms);

    const ContinuedFraction
      lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));

    //extrapolate the curve start
    if (lower_cf.is_valid()){
      double sample_size = static_cast<double>(sample);
      while(sample_size < max_extrapolation){
        double t = (sample_size - sample_vals_sum)/sample_vals_sum;
        assert(t >= 0.0);
        yield_vector.push_back(initial_distinct + t*lower_cf(t));
        sample_size += bin_step_size;
      }

      // SANITY CHECK
      if (check_yield_estimates(yield_vector)) {
        bootstrap_estimates.push_back(yield_vector);
        if (VERBOSE) cerr << '.';
      }
      else if (VERBOSE){
        cerr << "_";
      }
    }
    else if (VERBOSE){
      cerr << "_";
    }

  }
  if (VERBOSE)
    cerr << endl;
  if (bootstrap_estimates.size() < bootstraps)
    throw SMITHLABException("too many defects in the approximation, consider running in defect mode");
}

static bool
extrap_single_estimate(const bool VERBOSE, const bool DEFECTS,
		       vector<double> &hist,
                       size_t max_terms, const int diagonal,
                       const double step_size, 
                       const double max_extrapolation,
                       vector<double> &yield_estimate) {

  yield_estimate.clear();
  double vals_sum = 0.0;
  for(size_t i = 0; i < hist.size(); i++)
    vals_sum += i*hist[i];
  const double initial_distinct 
    = accumulate(hist.begin(), hist.end(), 0.0);

  // interpolate complexity curve by random sampling w/out replacement
  size_t upper_limit = static_cast<size_t>(vals_sum);
  size_t step = static_cast<size_t>(step_size);
  size_t sample = step;
  while (sample < upper_limit){
    yield_estimate.push_back(
		interpolate_distinct(hist, upper_limit, 
		                      static_cast<size_t>(initial_distinct), sample));
    sample += step;
  }

  // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
  size_t counts_before_first_zero = 1;
  while (counts_before_first_zero < hist.size() &&
         hist[counts_before_first_zero] > 0)
    ++counts_before_first_zero;


  // Ensure we are not using a zero term
  max_terms = std::min(max_terms, counts_before_first_zero - 1);

  // refit curve for lower bound (degree of approx is 1 less than
  // max_terms)
  max_terms = max_terms - (max_terms % 2 == 1);

  const ContinuedFractionApproximation
    lower_cfa(diagonal, max_terms);

  const ContinuedFraction
    lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));

  // extrapolate curve
  if (lower_cf.is_valid() || DEFECTS){
    double sample_size = static_cast<double>(sample);
    while(sample_size < max_extrapolation){
      const double one_minus_fold_extrap 
        = (sample_size - vals_sum)/vals_sum;
      assert(one_minus_fold_extrap >= 0.0);
      double tmp = one_minus_fold_extrap*lower_cf(one_minus_fold_extrap);
      yield_estimate.push_back(initial_distinct + tmp);
      sample_size += step_size;
    }
  }
  else{
    // FAIL!
    // lower_cf unacceptable, need to bootstrap to obtain estimates
    return false;
  }

  if (VERBOSE) {
    if(lower_cf.offset_coeffs.size() > 0){
      cerr << "CF_OFFSET_COEFF_ESTIMATES" << endl;
      copy(lower_cf.offset_coeffs.begin(), lower_cf.offset_coeffs.end(),
           std::ostream_iterator<double>(cerr, "\n"));
    }
    if(lower_cf.cf_coeffs.size() > 0){
      cerr << "CF_COEFF_ESTIMATES" << endl;
      copy(lower_cf.cf_coeffs.begin(), lower_cf.cf_coeffs.end(),
           std::ostream_iterator<double>(cerr, "\n"));
    }
  }

  // SUCCESS!!
  return true;
}

static double
GoodToulmin2xExtrap(const vector<double> &counts_hist){
  double two_fold_extrap = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    two_fold_extrap += pow(-1.0, i + 1)*counts_hist[i];

  return two_fold_extrap;
}


static void
write_predicted_complexity_curve(const string outfile,
                                 const double c_level, const double step_size,
                                 const vector<double> &yield_estimates,
                                 const vector<double> &yield_lower_ci_lognormal,
                                 const vector<double> &yield_upper_ci_lognormal) {
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
        << yield_lower_ci_lognormal[i] << '\t'
        << yield_upper_ci_lognormal[i] << endl;
}

static void
write_predicted_coverage_curve(const string outfile,
                               const double c_level,
                               const double base_step_size,
                               const size_t bin_size,
                               const vector<double> &coverage_estimates,
                               const vector<double> &coverage_lower_ci_lognormal,
                               const vector<double> &coverage_upper_ci_lognormal) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "TOTAL_BASES\tEXPECTED_COVERED_BASES\t"
      << "LOWER_" << 100*c_level << "%CI\t"
      << "UPPER_" << 100*c_level << "%CI" << endl;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < coverage_estimates.size(); ++i)
    out << (i + 1)*base_step_size << '\t'
        << coverage_estimates[i]*bin_size << '\t'
        << coverage_lower_ci_lognormal[i]*bin_size << '\t'
        << coverage_upper_ci_lognormal[i]*bin_size << endl;
}


static int
lc_extrap(const int argc, const char **argv) {
  
  try {
    const size_t MIN_REQUIRED_COUNTS = 4;
    
    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 100;
    double max_extrapolation = 1.0e10;
    double step_size = 1e6;
    size_t bootstraps = 100;
    int diagonal = 0;
    double c_level = 0.95;
    double dupl_level = 0.5;
      
    /* FLAGS */
    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool SINGLE_ESTIMATE = false;
    bool DEFECTS = false;
      
#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif
      
    /********** GET COMMAND LINE ARGUMENTS  FOR LC EXTRAP ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation "
                      "(default: " + toa(max_extrapolation) + ")",
                      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size in extrapolations "
                      "(default: " + toa(step_size) + ")",
                      false, step_size);
    opt_parse.add_opt("bootstraps",'n',"number of bootstraps "
                      "(default: " + toa(bootstraps) + "), ",
                      false, bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals "
                      "(default: " + toa(c_level) + ")", false, c_level);
    opt_parse.add_opt("dupl_level",'d', "fraction of duplicate to predict "
                      "(default: " + toa(dupl_level) + ")",
                      false, dupl_level);
    opt_parse.add_opt("terms",'x',"maximum number of terms", false,
                      orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
#ifdef HAVE_SAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format",
                      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads (default: "
                      + toa(MAX_SEGMENT_LENGTH) + ")",
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
    opt_parse.add_opt("quick",'Q',
                      "quick mode, estimate yield without bootstrapping for confidence intervals",
                      false, SINGLE_ESTIMATE);
    opt_parse.add_opt("defects", 'D', 
		      "defects mode to extrapolate without testing for defects",
		      false, DEFECTS);

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


    vector<double> counts_hist;
    size_t n_reads = 0;

    // LOAD VALUES
    if(HIST_INPUT){
      if(VERBOSE)
        cerr << "HIST_INPUT" << endl;
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if(VALS_INPUT){
      if(VERBOSE)
        cerr << "VALS_INPUT" << endl;
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_SAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BAM_INPUT" << endl;
      const size_t MAX_READS_TO_HOLD = 5000000;
      size_t n_paired = 0;
      size_t n_mates = 0;
      n_reads = load_counts_BAM_pe(VERBOSE, input_file_name, 
                                   MAX_SEGMENT_LENGTH, 
                                   MAX_READS_TO_HOLD, n_paired, 
                                   n_mates, counts_hist);
      if(VERBOSE){
        cerr << "MERGED PAIRED END READS = " << n_paired << endl;
        cerr << "MATES PROCESSED = " << n_mates << endl;
      }
    }
    else if(BAM_FORMAT_INPUT){
      if(VERBOSE)
        cerr << "BAM_INPUT" << endl;
      n_reads = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if(PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BED_INPUT" << endl;
      n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else{ // default is single end bed file
      if(VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }

    const size_t max_observed_count = counts_hist.size() - 1;
    const double distinct_reads = accumulate(counts_hist.begin(),
                                             counts_hist.end(), 0.0);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < counts_hist.size() &&
           counts_hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;

    orig_max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);


    const size_t distinct_counts =
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
                                        bind2nd(std::greater<double>(), 0.0)));
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
    if(two_fold_extrap < 0.0)
      throw SMITHLABException("Library expected to saturate in doubling of "
                              "size, unable to extrapolate");


    size_t total_reads = 0;
    for(size_t i = 0; i < counts_hist.size(); i++){
      total_reads += i*counts_hist[i];
    }
    //assert(total_reads == n_reads);

    // catch if all reads are distinct
    if (orig_max_terms < MIN_REQUIRED_COUNTS)
      throw SMITHLABException("max count before zero is les than min required "
                              "count (4), sample not sufficiently deep or "
                              "duplicates removed");

    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // ESTIMATE COMPLEXITY CURVE

    if(VERBOSE)
      cerr << "[ESTIMATING YIELD CURVE]" << endl;
    vector<double> yield_estimates;


    if(SINGLE_ESTIMATE || DEFECTS){
      bool SINGLE_ESTIMATE_SUCCESS =
        extrap_single_estimate(VERBOSE, DEFECTS, counts_hist, orig_max_terms, 
                               diagonal, step_size, max_extrapolation, 
                               yield_estimates);
      // IF FAILURE, EXIT
      if(!SINGLE_ESTIMATE_SUCCESS)
        throw SMITHLABException("SINGLE ESTIMATE FAILED, NEED TO RUN "
                                "FULL MODE FOR ESTIMATES");

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS\tEXPECTED_DISTINCT" << endl;

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << 0 << '\t' << 0 << endl;
      for (size_t i = 0; i < yield_estimates.size(); ++i)
        out << (i + 1)*step_size << '\t'
            << yield_estimates[i] << endl;

    }
    else{
      if (VERBOSE)
        cerr << "[BOOTSTRAPPING HISTOGRAM]" << endl;

      const size_t max_iter = 10*bootstraps;

      vector<vector <double> > bootstrap_estimates;
      extrap_bootstrap(VERBOSE, counts_hist, bootstraps, orig_max_terms,
                       diagonal, step_size, max_extrapolation, max_iter,
                       bootstrap_estimates);


      /////////////////////////////////////////////////////////////////////
      if (VERBOSE)
        cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;

      // yield ci
      vector<double> yield_upper_ci_lognormal, yield_lower_ci_lognormal;

      vector_median_and_ci(bootstrap_estimates, c_level, yield_estimates,
                           yield_lower_ci_lognormal, yield_upper_ci_lognormal);

      /////////////////////////////////////////////////////////////////////
      if (VERBOSE)
        cerr << "[WRITING OUTPUT]" << endl;

      write_predicted_complexity_curve(outfile, c_level, step_size,
                                       yield_estimates, yield_lower_ci_lognormal,
                                       yield_upper_ci_lognormal);
    }
  }
  catch (SMITHLABException &e) {
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
    double max_extrapolation = 1.0e12;
    size_t bootstraps = 100;
    bool DEFECTS = false;

    bool NO_SEQUENCE = false;
    double c_level = 0.95;

    // ********* GET COMMAND LINE ARGUMENTS  FOR GC EXTRAP **********
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "<sorted-mapped-read-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("max_width", 'w', "max fragment length, "
                      "set equal to read length for single end reads",
                      false, max_width);
    opt_parse.add_opt("bin_size", 'b', "bin size "
                      "(default: " + toa(bin_size) + ")",
                      false, bin_size);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation in base pairs"
                      "(default: " + toa(max_extrapolation) + ")",
                      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size in bases between extrapolations "
                      "(default: " + toa(base_step_size) + ")",
                      false, base_step_size);
    opt_parse.add_opt("bootstraps",'n',"number of bootstraps "
                      "(default: " + toa(bootstraps) + "), ",
                      false, bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals "
                      "(default: " + toa(c_level) + ")", false, c_level);
    opt_parse.add_opt("terms",'x',"maximum number of terms",
                      false, orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.add_opt("bed", 'D',
                      "input is in bed format without sequence information",
                      false, NO_SEQUENCE);
    opt_parse.add_opt("quick",'Q',
                      "quick mode: run gc_extrap without "
                      "bootstrapping for confidence intervals",
                      false, SINGLE_ESTIMATE);
    opt_parse.add_opt("defects", 'D', 
		      "defects mode to extrapolate without testing for defects",
		      false, DEFECTS);


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
    if(VERBOSE)
      cerr << "LOADING READS" << endl;

    if(NO_SEQUENCE){
      if(VERBOSE)
        cerr << "BED FORMAT" << endl;
      n_reads = load_coverage_counts_GR(input_file_name, bin_size,
                                        max_width, coverage_hist);
    }
    else{
      if(VERBOSE)
        cerr << "MAPPED READ FORMAT" << endl;
      n_reads = load_coverage_counts_MR(VERBOSE, input_file_name, bin_size,
                                        max_width, coverage_hist);
    }

    double total_bins = 0.0;
    for(size_t i = 0; i < coverage_hist.size(); i++)
      total_bins += coverage_hist[i]*i;
    const double distinct_bins = 
      accumulate(coverage_hist.begin(), coverage_hist.end(), 0.0);
    
    const double avg_bins_per_read = total_bins/n_reads;
    double bin_step_size = base_step_size/bin_size;

    const size_t max_observed_count = coverage_hist.size() - 1;

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < coverage_hist.size() &&
           coverage_hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;

    orig_max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);

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
      throw SMITHLABException("max count before zero is les than min required "
                              "count (4), sample not sufficiently deep or "
                              "duplicates removed");

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(coverage_hist);
    if(two_fold_extrap < 0.0)
      throw SMITHLABException("Library expected to saturate in doubling of "
                              "experiment size, unable to extrapolate");


    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // ESTIMATE COMPLEXITY CURVE

    if(VERBOSE)
      cerr << "[ESTIMATING COVERAGE CURVE]" << endl;
    vector<double> coverage_estimates;


    if (SINGLE_ESTIMATE) {
      
      bool SINGLE_ESTIMATE_SUCCESS =
        extrap_single_estimate(VERBOSE, DEFECTS, coverage_hist, orig_max_terms, diagonal,
                               bin_step_size, max_extrapolation/bin_size,
                               coverage_estimates);
      
      
      // IF FAILURE, EXIT
      if (!SINGLE_ESTIMATE_SUCCESS)
        throw SMITHLABException("SINGLE ESTIMATE FAILED, NEED TO RUN IN "
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
      
      const size_t max_iter = 10*bootstraps;
      
      vector<vector <double> > bootstrap_estimates;
      extrap_bootstrap(VERBOSE, coverage_hist, bootstraps, orig_max_terms,
                       diagonal, bin_step_size, max_extrapolation/bin_size,
                       max_iter, bootstrap_estimates);
      
      
      /////////////////////////////////////////////////////////////////////
      if (VERBOSE)
        cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      
      vector<double> coverage_upper_ci_lognormal, coverage_lower_ci_lognormal;
      vector_median_and_ci(bootstrap_estimates, c_level, 
                           coverage_estimates, coverage_lower_ci_lognormal, 
                           coverage_upper_ci_lognormal);
      
      /////////////////////////////////////////////////////////////////////
      if (VERBOSE)
        cerr << "[WRITING OUTPUT]" << endl;
      write_predicted_coverage_curve(outfile, c_level, base_step_size,
                                     bin_size, coverage_estimates,
                                     coverage_lower_ci_lognormal,
                                     coverage_upper_ci_lognormal);
    }
  }
  catch (SMITHLABException &e) {
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
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////  C_CURVE BELOW HERE
/////

static int
c_curve(const int argc, const char **argv) {

  try {  

    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;

    string outfile;

    size_t upper_limit = 0;
    double step_size = 1e6;
  
#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("step",'s',"step size in extrapolations "
                      "(default: " + toa(step_size) + ")",
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
#ifdef HAVE_SAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format",
                      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads (default: "
                      + toa(MAX_SEGMENT_LENGTH) + ")",
                      false, MAX_SEGMENT_LENGTH);
#endif
  
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
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); // use default type
    srand(time(0) + getpid()); //give the random fxn a new seed
    gsl_rng_set(rng, rand()); //initialize random number generator with the seed

    vector<double> counts_hist;
    size_t n_reads = 0;

    // LOAD VALUES
    if(HIST_INPUT){
      if(VERBOSE)
        cerr << "INPUT_HIST" << endl;
      n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if (VALS_INPUT) {
      if (VERBOSE)
        cerr << "VALS_INPUT" << endl;
      n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_SAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END){
      if(VERBOSE)
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
    const double distinct_reads = accumulate(counts_hist.begin(),
                                             counts_hist.end(), 0.0);
  
    size_t total_reads = 0;
    for(size_t i = 0; i < counts_hist.size(); i++)
      total_reads += i*counts_hist[i];
  
    const size_t distinct_counts =
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
                                        bind2nd(std::greater<double>(), 0.0)));

    if (VERBOSE)
      cerr << "TOTAL READS     = " << n_reads << endl
           << "COUNTS_SUM      = " << total_reads << endl
           << "DISTINCT READS  = " << distinct_reads << endl
           << "DISTINCT COUNTS = " << distinct_counts << endl
           << "MAX COUNT       = " << max_observed_count << endl
           << "COUNTS OF 1     = " << counts_hist[1] << endl;
  
  
    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
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
  catch (SMITHLABException &e) {
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

  try{

    bool VERBOSE = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool VALS_INPUT = false;
    bool QUICK_MODE = false;

    string outfile;

#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 5000;
#endif

    size_t max_num_points = 10;
    double tolerance = 1e-20;
    size_t bootstraps = 100;
    double c_level = 0.95;
    size_t max_iter = 100;


    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("max_num_points",'p',"maximum number of points in quadrature "
                      "estimates (default: " + toa(max_num_points) + ")",
                      false, max_num_points);
    opt_parse.add_opt("tolerance", 't', "numerical tolerance "
                      "(default: " + toa(tolerance) + ")",
		      false, tolerance);
    opt_parse.add_opt("bootstraps", 'n', "number of bootstraps "
                      "(default: " + toa(bootstraps) + ")",
		      false, bootstraps);
    opt_parse.add_opt("clevel", 'c', "level for confidence intervals "
                      "(default: " + toa(c_level) + ")", false, c_level);
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
#ifdef HAVE_SAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format",
                      false, BAM_FORMAT_INPUT);
    opt_parse.add_opt("seg_len", 'l', "maximum segment length when merging "
                      "paired end bam reads (default: "
                      + toa(MAX_SEGMENT_LENGTH) + ")",
                      false, MAX_SEGMENT_LENGTH);
#endif
    opt_parse.add_opt("quick", 'q', "quick mode, estimates without bootstrapping",
		      false, QUICK_MODE);


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
    if(HIST_INPUT){
      if(VERBOSE)
        cerr << "HIST_INPUT" << endl;
      n_obs = load_histogram(input_file_name, counts_hist);
    }
    else if(VALS_INPUT){
      if(VERBOSE)
        cerr << "VALS_INPUT" << endl;
      n_obs = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_SAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BAM_INPUT" << endl;
      const size_t MAX_READS_TO_HOLD = 5000000;
      size_t n_paired = 0;
      size_t n_mates = 0;
      n_obs = load_counts_BAM_pe(VERBOSE, input_file_name, 
                                   MAX_SEGMENT_LENGTH, 
                                   MAX_READS_TO_HOLD, n_paired, 
                                   n_mates, counts_hist);
      if(VERBOSE){
        cerr << "MERGED PAIRED END READS = " << n_paired << endl;
        cerr << "MATES PROCESSED = " << n_mates << endl;
      }
    }
    else if(BAM_FORMAT_INPUT){
      if(VERBOSE)
        cerr << "BAM_INPUT" << endl;
      n_obs = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if(PAIRED_END){
      if(VERBOSE)
        cerr << "PAIRED_END_BED_INPUT" << endl;
      n_obs = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else{ // default is single end bed file
      if(VERBOSE)
        cerr << "BED_INPUT" << endl;
      n_obs = load_counts_BED_se(input_file_name, counts_hist);
    }

    const double distinct_obs = accumulate(counts_hist.begin(), 
					   counts_hist.end(), 0.0);


    vector<double> measure_moments;
    // mu_r = (r + 1)! n_{r+1} / n_1
    size_t indx = 1;
    while(counts_hist[indx] > 0  && indx <= counts_hist.size()){
      measure_moments.push_back(exp(gsl_sf_lnfact(indx)
				    + log(counts_hist[indx])
				    - log(counts_hist[1])));
      if(!std::isfinite(measure_moments.back())){
	measure_moments.pop_back();
	break;
      }
      indx++;
    }
    

    if (VERBOSE){
      cerr << "TOTAL OBSERVATIONS     = " << n_obs << endl
           << "DISTINCT OBSERVATIONS  = " << distinct_obs << endl
           << "MAX COUNT              = " << counts_hist.size() - 1 << endl;

      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << setprecision(16) << counts_hist[i] << endl;

      cerr << "OBSERVED MOMENTS" << endl;
      for(size_t i = 0; i < measure_moments.size(); i++)
	cerr << std::setprecision(16) << measure_moments[i] << endl;  
    }


    if(QUICK_MODE){
      if(measure_moments.size() < 2*max_num_points)
	max_num_points = static_cast<size_t>(floor(measure_moments.size()/2));
      else
	measure_moments.resize(2*max_num_points);
      size_t n_points = 0;
      n_points = ensure_pos_def_mom_seq(measure_moments, tolerance, VERBOSE);
      if(VERBOSE)
	cerr << "n_points = " << n_points << endl;    

      MomentSequence obs_mom_seq(measure_moments);
    
      if(VERBOSE){
	for(size_t k = 0; k < obs_mom_seq.alpha.size(); k++)
	  cerr << "alpha_" << k << '\t';
	cerr << endl;
	for(size_t k = 0; k < obs_mom_seq.alpha.size(); k++)
	  cerr << obs_mom_seq.alpha[k] << '\t';
	cerr << endl;

	for(size_t k = 0; k < obs_mom_seq.beta.size(); k++)
	  cerr << "beta_" << k << '\t';
	cerr << endl;
	for(size_t k = 0; k < obs_mom_seq.beta.size(); k++)
	  cerr << obs_mom_seq.beta[k] << '\t';
	cerr << endl;
      }
    
      vector<double> points, weights;
      obs_mom_seq.Lower_quadrature_rules(VERBOSE, n_points, tolerance,
					 max_iter, points, weights);

      const double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
      if(weights_sum != 1.0){
	for(size_t i = 0; i < weights.size(); i++)
	  weights[i] = weights[i]/weights_sum;
      }

      if(VERBOSE){
	cerr << "points = " << endl;
	for(size_t i = 0; i < points.size(); i++)
	  cerr << points[i] << '\t';
	cerr << endl;

	cerr << "weights = " << endl;
	for(size_t i = 0; i < weights.size(); i++)
	  cerr << weights[i] << '\t';
	cerr << endl;
      }

      double estimated_unobs = 0.0;
    
      for(size_t i = 0; i < weights.size(); i++)
	estimated_unobs += counts_hist[1]*weights[i]/points[i];

      if(estimated_unobs > 0.0)
	estimated_unobs += distinct_obs;
      else{
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
   else{
      vector<double> quad_estimates;

  //setup rng
      srand(time(0) + getpid());
      gsl_rng_env_setup();
      gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(rng, rand());

      // hist may be sparse, to speed up bootstrapping
      // sample only from positive entries
      vector<size_t> counts_hist_distinct_counts;
      vector<double> distinct_counts_hist;
      for (size_t i = 0; i < counts_hist.size(); i++){
	if (counts_hist[i] > 0) {
	  counts_hist_distinct_counts.push_back(i);
	  distinct_counts_hist.push_back(counts_hist[i]);
	}
      }

      for(size_t iter = 0; 
	  iter < max_iter && quad_estimates.size() < bootstraps; 
	  ++iter){

	vector<double> sample_hist;
	resample_hist(rng, counts_hist_distinct_counts, 
		      distinct_counts_hist, sample_hist);

	const double sampled_distinct = accumulate(sample_hist.begin(), sample_hist.end(), 0.0);
	// initialize moments, 0th moment is 1
	vector<double> bootstrap_moments(1, 1.0);
	// moments[r] = (r + 1)! n_{r+1} / n_1
	for(size_t i = 0; i < 2*max_num_points + 1; i++)
	  bootstrap_moments.push_back(exp(gsl_sf_lnfact(i + 2) 
					  + log(sample_hist[i + 2])
					  - log(sample_hist[1])) );

	size_t n_points = 0;
	n_points = ensure_pos_def_mom_seq(bootstrap_moments, tolerance, false);


	MomentSequence bootstrap_mom_seq(bootstrap_moments);

   	vector<double> points, weights;
	bootstrap_mom_seq.Lower_quadrature_rules(VERBOSE, n_points, tolerance,
						 max_iter, points, weights);

	const double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
	if(weights_sum != 1.0){
	  for(size_t i = 0; i < weights.size(); i++)
	    weights[i] = weights[i]/weights_sum;
	}

	double estimated_unobs = 0.0;
    
	for(size_t i = 0; i < weights.size(); i++)
	  estimated_unobs += counts_hist[1]*weights[i]/points[i];

	if(estimated_unobs > 0.0)
	  estimated_unobs += sampled_distinct;
	else{
	  estimated_unobs = sampled_distinct;
	  n_points = 0;
	}

	quad_estimates.push_back(estimated_unobs);
      }

      double log_mean_estimate, lower_log_ci, upper_log_ci;

      log_mean(VERBOSE, quad_estimates, c_level, log_mean_estimate, 
	       lower_log_ci, upper_log_ci);

     std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(1);

      out << "log_mean_estimated_unobs" << '\t'
	  << "log_lower_ci" << '\t'
	  << "log_upper_ci" << endl;
      out << log_mean_estimate << '\t'
	  << lower_log_ci << '\t'
	  << upper_log_ci << endl;


   }


  }
  catch (SMITHLABException &e) {
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
  
  static const string 
    USAGE_MESSAGE("preseq:  a program for analyzing library complexity\n"
                  "Version: " + toa(PRESEQ_VERSION) + "\n\n"
                  "Usage:   preseq <command> [OPTIONS]\n\n"
                  "<command>: c_curve    generate complexity curve for a library\n"
                  "           lc_extrap  predict the yield for future experiments\n"
                  "           gc_extrap  predict genome coverage low input\n"
                  "                      sequencing experiments\n"
		  "           bound_pop  lower bound on population size\n"
                  );
  
  if (argc < 2)
    cerr << USAGE_MESSAGE << endl;

  else if (strcmp(argv[1], "lc_extrap") == 0) {
    
    return lc_extrap(argc, argv);
    
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
  else {
    cerr << "unrecognized command: " << argv[1] << endl
         << USAGE_MESSAGE << endl;
    return EXIT_SUCCESS;
    
  }
}
