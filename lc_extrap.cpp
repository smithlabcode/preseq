/*    lc_extrap:
 *
 *    Copyright (C) 2012 University of Southern California and
 *                       Andrew D. Smith and Timothy Daley
 *
 *    Authors: Andrew D. Smith and Timothy Daley
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
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "pade_approximant.hpp"
#include "continued_fraction.hpp"
#include "library_size_estimates.hpp"

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <RNG.hpp>

#include <smithlab_os.hpp>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <fstream>
#include <numeric>
#include <vector>

using std::string;
using std::vector;
using std::count_if;
using std::endl;
using std::cerr;
using std::max;
using std::ostream;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;

#ifdef HAVE_BAMTOOLS
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;


static GenomicRegion
BamToGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
		   const BamAlignment &ba) {
  const unordered_map<size_t, string>::const_iterator
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));

  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = start + ba.Length;
  return GenomicRegion(chrom, start, end, "X", 0, ba.strand);
}


static void
load_values_BAM(const string &infile, vector<double> &values) {

  BamReader reader;
  reader.Open(infile);

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;

  size_t n_reads = 1;
  values.push_back(1.0);

  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    const GenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam));
    if (r < prev)
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start() ||
	r.get_strand() != prev.get_strand())
      values.push_back(1.0);
    else values.back()++;
    ++n_reads;
    prev.swap(r);
  }
  reader.Close();
  return n_reads;
}
#endif


static size_t
load_values(const string input_file_name, vector<double> &values) {

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw "problem opening file: " + input_file_name;

  GenomicRegion r, prev;
  if (!(in >> prev))
    throw "problem reading from: " + input_file_name;

  size_t n_reads = 1;
  values.push_back(1.0);
  while (in >> r) {
    if (r < prev)
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start() ||
	r.get_strand() != prev.get_strand())
      values.push_back(1.0);
    else values.back()++;
    ++n_reads;
    prev.swap(r);
  }
  return n_reads;
}


void
resample_hist(const vector<double> &vals_hist,
	      const gsl_rng *rng,
	      const double total_sampled_reads,
	      double expected_sample_size,
	      vector<double> &sample_hist) {
  sample_hist.resize(vals_hist.size());
  double remaining = total_sampled_reads;
  const double vals_mean = total_sampled_reads/expected_sample_size;
 
  expected_sample_size /= 2.0;
  vector<unsigned int> tmp(vals_hist.size());
  gsl_ran_multinomial(rng, vals_hist.size(), expected_sample_size,
		      &vals_hist.front(), &tmp.front());
  double inc = 0.0;
  for (size_t i = 0; i < tmp.size(); ++i)
    inc += i*tmp[i];

  vector<double> double_tmp(tmp.size());
  for (size_t i = 0; i < tmp.size(); i++)
    double_tmp[i] = static_cast<double>(tmp[i]);
  // add the new stuff
  transform(double_tmp.begin(), double_tmp.end(), sample_hist.begin(), 
	    sample_hist.begin(), std::plus<double>());
  // update the amount we still need to get
  remaining -= inc;

  while (remaining > 0 && inc > 0) {

    // get a new sample
    expected_sample_size = remaining/vals_mean;
    expected_sample_size /= 2.0;
    tmp.clear();
    tmp.resize(vals_hist.size());
    gsl_ran_multinomial(rng, vals_hist.size(), expected_sample_size,
			&vals_hist.front(), &tmp.front());
    // see how much we got
    inc = 0.0;
    for (size_t i = 0; i < tmp.size(); ++i)
      inc += i*tmp[i];

    // only add to histogram if sampled reads < remaining reads
    if (inc <= remaining) {
      double_tmp.clear();
      double_tmp.resize(tmp.size());
      for (size_t i = 0; i < tmp.size(); i++)
	double_tmp[i] = static_cast<double>(tmp[i]);
      // add the new stuff
      transform(double_tmp.begin(), double_tmp.end(), sample_hist.begin(),
		sample_hist.begin(), std::plus<double>());
      // update the amount we still need to get
      remaining -= inc;
    }

  }

  // just right
  assert(remaining >= 0);
  sample_hist[static_cast<size_t>(remaining)]++;

}


static bool
check_yield_estimates(const vector<double> &estimates) {
  
  if (estimates.empty()) 
    return false;

  // make sure that the estimate is increasing in the time_step and is
  // below the initial distinct per step_size
  if (!finite(accumulate(estimates.begin(), estimates.end(), 0.0)))
    return false;
  
  for (size_t i = 1; i < estimates.size(); ++i)
    if ((estimates[i] < estimates[i - 1]) ||
	(i >= 2 && (estimates[i] - estimates[i - 1] >
		    estimates[i - 1] - estimates[i - 2])) ||
	(estimates[i] < 0.0))
      return false;
  
  return true;
}

static inline bool
check_saturation_estimates(const vector<double> estimates) {
  if (estimates.empty())
    return false;

  // make sure estimates are decreasing and
  // between 0 & 1
  if (estimates[0] >= 1.0 || estimates[0] < 0.0)
    return false;

  for(size_t i = 1; i < estimates.size(); i++)
    if (estimates[i] > estimates[i-1] ||
	estimates[i] >= 1.0 ||
	estimates[i] < 0.0) 
      return false;
  
  return true;
}

void
estimates_bootstrap(const bool VERBOSE, const vector<double> &orig_values, 
		    const size_t bootstraps, const size_t orig_max_terms, 
		    const int diagonal, const double step_size, 
		    const double max_extrapolation, 
		    const double max_val, const double val_step, 
		    vector<double> &lower_bound_size,
		    vector<double> &upper_bound_size,
		    vector< vector<double> > &saturation_estimates,
		    vector< vector<double> > &yield_estimates) {
  // clear returning vectors
  upper_bound_size.clear();
  lower_bound_size.clear();
  yield_estimates.clear();
  saturation_estimates.clear();
  
  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand()); 

  const size_t max_observed_count = 
    static_cast<size_t>(*std::max_element(orig_values.begin(), 
					  orig_values.end()));
    
  vector<double> orig_hist(max_observed_count + 1, 0.0);
  for (size_t i = 0; i < orig_values.size(); ++i)
    ++orig_hist[static_cast<size_t>(orig_values[i])];


  const double vals_sum = accumulate(orig_values.begin(), orig_values.end(), 0.0);

  for (size_t iter = 0; 
       iter < 2*bootstraps && yield_estimates.size() < bootstraps; ++iter) {
    
    vector<double> hist;
    resample_hist(orig_hist, rng, vals_sum,
		  static_cast<double>(orig_values.size()), hist);
    
    //resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();
    
    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < hist.size() &&
	   hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
    
    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 0);
    
    //refit curve for lower bound
    const ContinuedFractionApproximation 
      lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
    
    const ContinuedFraction 
      lower_cf(lower_cfa.optimal_cont_frac_yield(hist));
    
    vector<double> yield_vector;
    vector<double> saturation_vector;
    if (lower_cf.is_valid()) {
      lower_cf.extrapolate_distinct(hist, max_val, 
				    val_step, yield_vector);
      lower_cf.extrapolate_yield_deriv(hist, vals_sum, max_val,
				       val_step, saturation_vector);
    }
    // SANITY CHECK
    const bool ACCEPT_ESTIMATES = 
      check_yield_estimates(yield_vector) && 
      check_saturation_estimates(saturation_vector);
    if (ACCEPT_ESTIMATES) {
      const double distinct = accumulate(hist.begin(), hist.end(), 0.0);
      const double upper_bound =
	upperbound_librarysize(false, hist, lower_cf.return_degree()) + distinct;
      const double lower_bound = 
	lower_cfa.lowerbound_librarysize(hist, upper_bound);
      
      if (finite(lower_bound) && lower_bound > 0 && finite(upper_bound)) {
	
	yield_estimates.push_back(yield_vector);
	saturation_estimates.push_back(saturation_vector);
	upper_bound_size.push_back(upper_bound);
	lower_bound_size.push_back(lower_bound);
	if (VERBOSE) cerr << '.';
      }
      else if (VERBOSE) cerr << '_';
    }
    else if (VERBOSE) cerr << '_';
    if (iter == 2*bootstraps - 1)
      throw SMITHLABException("too many iterations, poor sample");
  }
  if (VERBOSE)
    cerr << endl;
}

void
single_estimates(bool VERBOSE, const vector<double> &hist,
		 size_t max_terms, const int diagonal, const double vals_sum,
		 const double step_size, const double max_extrapolation, 
		 const double max_val, const double val_step,
		 double &upper_bound, double &lower_bound,
		 vector<double> &saturation_estimates,
		 vector<double> &yield_estimates) {

  saturation_estimates.clear();
  yield_estimates.clear();
  
  // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
  size_t counts_before_first_zero = 1;
  while (counts_before_first_zero < hist.size() &&
	 hist[counts_before_first_zero] > 0)
    ++counts_before_first_zero;
    
  max_terms = std::min(max_terms, counts_before_first_zero - 1);
  // refit curve for lower bound (degree of approx is 1 less than
  // max_terms)
  max_terms = max_terms - (max_terms % 2 == 0);

  //refit curve for lower bound
  const ContinuedFractionApproximation 
    lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
    
  const ContinuedFraction 
    lower_cf(lower_cfa.optimal_cont_frac_yield(hist));
    
  if (lower_cf.is_valid()) {
    lower_cf.extrapolate_distinct(hist, max_val, 
				  val_step, yield_estimates);
    lower_cf.extrapolate_yield_deriv(hist, vals_sum, max_val,
				     val_step, saturation_estimates);
  }
  if (VERBOSE) {
    cerr << "ESTIMATES_CONTINUED_FRACTION_OFFSET_COEFFS" << endl;
    vector<double> off_coeffs(lower_cf.offset_coeffs);
    for (size_t j = 0; j < off_coeffs.size(); j++)
      cerr << off_coeffs[j] << endl;
    cerr << "ESTIMATES_CONTINUED_FRACTION_COEFFS" << endl;
    vector<double> cf_coefs(lower_cf.cf_coeffs);
    for (size_t j = 0; j < cf_coefs.size(); j++)
      cerr << cf_coefs[j] << endl;
  }
    
  // SANITY CHECK
  if (check_yield_estimates(yield_estimates) && 
      check_saturation_estimates(saturation_estimates)) {
      
    const double distinct = accumulate(hist.begin(), hist.end(), 0.0);
    upper_bound =
      upperbound_librarysize(VERBOSE, hist, lower_cf.return_degree()) + distinct;
    
    ContinuedFraction lb_opt;
    lower_bound = 
      lower_cfa.lowerbound_librarysize(hist, upper_bound, lb_opt);
    
    if (VERBOSE) {
      cerr << "LOWER_BOUND_CONTINUED_FRACTION_OFFSET_COEFFS" << endl;
      vector<double> off_coeffs(lb_opt.offset_coeffs);
      for(size_t j = 0; j < off_coeffs.size(); j++)
	cerr << off_coeffs[j] << endl;
      cerr << "LOWER_BOUND_CONTINUED_FRACTION_COEFF" << endl;
      vector<double> cf_coefs(lb_opt.cf_coeffs);
      for(size_t j = 0; j < cf_coefs.size(); j++)
	cerr << cf_coefs[j] << endl;
    }
    
    // error if bounds are unacceptable
    if (!finite(upper_bound))
      cerr << "UPPER_BOUND INFINITE" << endl;
    if (!finite(lower_bound) || lower_bound <= 0.0)
      cerr << "LOWER_BOUND UNACCEPTABLE" << endl;
  }
  else
    cerr << "ESTIMATES UNSTABLE, MORE DATA REQUIRED" << endl;


}



static double
compute_var(const vector<double> &estimates) {
  const double sample_size = estimates.size();
  const double mean = 
    accumulate(estimates.begin(), estimates.end(), 0.0)/sample_size;
  double variance = 0.0;
  for (size_t i = 0; i < estimates.size(); i++)
    variance += (estimates[i] - mean)*(estimates[i] - mean)/sample_size;
  return variance; 
}


static inline double
alpha_log_confint_multiplier(const double estimate,
			     const double initial_distinct,
			     const double variance,
			     const double alpha) {
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  return exp(inv_norm_alpha*
	     sqrt(log(1.0 + variance/pow(estimate - 
					 initial_distinct, 2))));
}

static void
return_median_and_alphaCI(const vector<vector<double> > &estimates,
			  const double alpha, const double initial_distinct,
			  vector<double> &median_estimates,
			  vector<double> &lower_alphaCI,
			  vector<double> &upper_alphaCI) {
  
  for (size_t i = 0; i < estimates[0].size(); i++) {
    // estimates is in wrong order, work locally on const val
    vector<double> estimates_row(estimates.size(), 0.0);
    for (size_t k = 0; k < estimates_row.size(); ++k)
      estimates_row[k] = estimates[k][i];
      
    // sort to get confidence interval
    sort(estimates_row.begin(), estimates_row.end());
    median_estimates.push_back(estimates_row[estimates.size()/2 - 1]);
    const double variance = compute_var(estimates_row);
    const double confint_multiplier = 
      alpha_log_confint_multiplier(median_estimates.back(), initial_distinct,
				   variance, alpha);
    upper_alphaCI.push_back(initial_distinct + 
			    (median_estimates.back() - 
			     initial_distinct)*confint_multiplier);
    lower_alphaCI.push_back(initial_distinct + 
			    (median_estimates.back() - 
			     initial_distinct)/confint_multiplier);
  }

}


static void
write_bootstrap_size_info(const string size_outfile,
			  const double c_level,
			  const double initial_distinct,
			  const vector<double> &lower_libsize,
			  const vector<double> &upper_libsize) {
  
  std::ofstream st_of;
  if (!size_outfile.empty()) 
    st_of.open(size_outfile.c_str());
  ostream size_out(size_outfile.empty() ? cerr.rdbuf() : st_of.rdbuf());
      
  const double cf_upper_size_var = compute_var(upper_libsize);
  const double cf_upper_size_mean = 
    accumulate(upper_libsize.begin(), 
	       upper_libsize.end(), 0.0)/upper_libsize.size();
      
  const double cf_upper_size_alpha_multiplier = 
    alpha_log_confint_multiplier(cf_upper_size_mean, initial_distinct, 
				 cf_upper_size_var, 1.0 - c_level);
      
  size_out << "LIBRARY_SIZE_UPPERBOUND_MEAN\t" 
	   << setprecision(0) << std::fixed 
	   << cf_upper_size_mean << endl;
      
  const double lib_size_ub_lower = initial_distinct + 
    (cf_upper_size_mean - initial_distinct)/cf_upper_size_alpha_multiplier;
  size_out << "LIBRARY_SIZE_UPPERBOUND_LOWER "
	   << setprecision(0) << 100*c_level << "% CI\t"
	   << fixed << lib_size_ub_lower << endl;

  const double lib_size_ub_upper = initial_distinct + 
    (cf_upper_size_mean - initial_distinct)*cf_upper_size_alpha_multiplier;
  size_out << "LIBRARY_SIZE_UPPERBOUND_UPPER "
	   << setprecision(0) << 100*c_level << "% CI\t"
	   << fixed << lib_size_ub_upper << endl;
      
  // NOW OUTPUT VALUES RELATED TO LIBRARY SIZE LOWER BOUND
  const double cf_lower_size_var = compute_var(lower_libsize);
  const double cf_lower_size_mean = 
    accumulate(lower_libsize.begin(), 
	       lower_libsize.end(), 0.0)/lower_libsize.size();
	
  const double cf_lower_size_alpha_multiplier = 
    alpha_log_confint_multiplier(cf_lower_size_mean, initial_distinct, 
				 cf_lower_size_var, 1.0 - c_level);
      
  size_out << "LIBRARY_SIZE_LOWERBOUND_MEAN\t" 
	   << fixed << cf_lower_size_mean << endl;
      
  const double lib_size_lb_lower = initial_distinct + 
    (cf_lower_size_mean - initial_distinct)/cf_lower_size_alpha_multiplier;
  size_out << "LIBRARY_SIZE_LOWERBOUND_LOWER " 
	   << setprecision(0) << 100*c_level << "% CI\t"
	   << fixed << lib_size_lb_lower << endl;
	
  const double lib_size_lb_upper = initial_distinct + 
    (cf_lower_size_mean - initial_distinct)*cf_lower_size_alpha_multiplier;
  size_out << "LIBRARY_SIZE_LOWERBOUND_UPPER " 
	   << setprecision(0) << 100*c_level << "% CI\t"
	   << fixed << lib_size_lb_upper << endl;
}



static void
write_predicted_curve(const string outfile, const double values_sum,
		      const double c_level, const double val_step,
		      const vector<double> &median_yield_estimates,
		      const vector<double> &yield_lower_alphaCI,
		      const vector<double> &yield_upper_alphaCI) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out << "#TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOWER_" << 100*c_level << "%CI\t"
      << "UPPER_" << 100*c_level << "%CI" << endl;
  
  double val = 0.0;
  for (size_t i = 0; i < median_yield_estimates.size(); ++i, val += val_step)
    out << fixed << setprecision(1) 
	<< (val + 1.0)*values_sum << '\t' 
	<< median_yield_estimates[i] << '\t'
	<< yield_lower_alphaCI[i] << '\t' 
	<< yield_upper_alphaCI[i] << endl;
}

int
main(const int argc, const char **argv) {

  try {

    const size_t MIN_REQUIRED_COUNTS = 5;

    /* FILES */
    string outfile;
    string size_outfile;
    string saturation_outfile;
    
    size_t orig_max_terms = 200;
    double max_extrapolation = 1.0e10;
    double step_size = 1e6;
    size_t bootstraps = 100;
    int diagonal = -1;
    double c_level = 0.95;
    
    /* FLAGS */
    bool VERBOSE = false;
    //  bool SMOOTH_HISTOGRAM = false;  
    
#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "",
			   "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
		      false , outfile);
    opt_parse.add_opt("LIBRARY_SIZE", 'L', "library size output file", 
		      false , size_outfile);
    opt_parse.add_opt("SATURATION", 'S', "saturation output file",
		      false, saturation_outfile);
    opt_parse.add_opt("extrapolation_length",'e',
		      "maximum extrapolation length "
		      "(default: " + toa(max_extrapolation) + ")",
                      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size in extrapolations "
		      "(default: " + toa(step_size) + ")", 
		      false, step_size);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps "
		      "(default: " + toa(bootstraps) + "), "
		      "set to 0 or 1 for single estimate",
		      false, bootstraps);
    opt_parse.add_opt("c_level", 'c', "level for confidence intervals "
		      "(default: " + toa(c_level) + ")", false, c_level);
    //  opt_parse.add_opt("terms",'t',"maximum number of terms", false, 
    //	     orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false, VERBOSE);
#ifdef HAVE_BAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format", 
		      false, BAM_FORMAT_INPUT);
#endif
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    
    vector<double> values;
#ifdef HAVE_BAMTOOLS
    if (BAM_FORMAT_INPUT)
      load_values_BAM(input_file_name, values);
    else
#endif
    load_values(input_file_name, values);

    const double initial_distinct = static_cast<double>(values.size());
    
    // JUST A SANITY CHECK
    const size_t values_sum = 
      static_cast<double>(accumulate(values.begin(), values.end(), 0.0));
    
    const double max_val = max_extrapolation/static_cast<double>(values_sum);
    const double val_step = step_size/static_cast<double>(values_sum);
    
    const size_t max_observed_count = 
      static_cast<size_t>(*std::max_element(values.begin(), values.end()));

    // catch if all reads are distinct
    if (max_observed_count < MIN_REQUIRED_COUNTS)
      throw SMITHLABException("sample appears too uniform");
    
    // BUILD THE HISTOGRAM
    vector<double> counts_hist(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_hist[static_cast<size_t>(values[i])];
    
    const size_t distinct_counts = 
      static_cast<size_t>(count_if(counts_hist.begin(), counts_hist.end(),
				   bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
      cerr << "TOTAL READS     = " << values_sum << endl
	   << "DISTINCT COUNTS = " << distinct_counts << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_hist[1] << endl
	   << "MAX TERMS       = " << orig_max_terms << endl;
    
    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    }
    
    // bootstrap the histogram, estimates are average
    if (bootstraps > 10) {
      if (VERBOSE) 
	cerr << "[BOOTSTRAP ESTIMATES]" << endl;
      
      vector<vector <double> > yield_estimates;
      vector< vector<double> > saturation_estimates;
      vector<double> lower_libsize, upper_libsize;
      estimates_bootstrap(VERBOSE, values,  bootstraps, orig_max_terms,
			  diagonal, step_size, max_extrapolation, 
			  max_val, val_step, lower_libsize, upper_libsize,
			  saturation_estimates, yield_estimates);
      
      if (VERBOSE)
	cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      
      vector<double> median_yield_estimates;
      vector<double> yield_upper_alphaCI, yield_lower_alphaCI;
      return_median_and_alphaCI(yield_estimates, 1.0 - c_level, 
				initial_distinct, median_yield_estimates, 
				yield_lower_alphaCI, yield_upper_alphaCI);

      vector<double> median_saturation_estimates;
      vector<double> saturation_upper_alphaCI, saturation_lower_alphaCI;
      return_median_and_alphaCI(saturation_estimates, 1.0 - c_level, 
				initial_distinct, median_saturation_estimates, 
				saturation_lower_alphaCI, 
				saturation_upper_alphaCI);
      
      if (VERBOSE) 
	cerr << "[WRITING OUTPUT]" << endl;
      
      write_predicted_curve(outfile, values_sum, c_level, val_step,
			    median_yield_estimates,
			    yield_lower_alphaCI, yield_upper_alphaCI);
      
      // IF VERBOSE OUTPUT IS REQUESTED, OR IF STATS ARE REQUESTED IN A
      // FILE, PRINT THEM!!
      if (VERBOSE || !size_outfile.empty())
	write_bootstrap_size_info(size_outfile, c_level, initial_distinct,
				  lower_libsize, upper_libsize);
      
      if (!saturation_outfile.empty()) {
	std::ofstream s_of;
	if (!saturation_outfile.empty())
	  s_of.open(saturation_outfile.c_str());
	ostream sat_out(s_of.rdbuf());
	sat_out << "#TOTAL_REANDS" << "\t"
		<< "SATURATION" << "\t"
		<< "LOWER_" << 100*c_level << "%CI" << "\t"
		<< "UPPER_" << 100*c_level << "%CI" << endl;
	double val = 0.0;
	for (size_t i = 0; i < median_saturation_estimates.size(); ++i, val += val_step)
	  sat_out << fixed << setprecision(1) 
		  << (val + 1.0)*values_sum << '\t' 
		  << median_saturation_estimates[i] << '\t'
		  << saturation_lower_alphaCI[i] << '\t' 
		  << saturation_upper_alphaCI[i] << endl;
      }
    }
    // if there are not enough bootstraps, do single estimate
    else{
      
      if (VERBOSE)
	cerr << "[ESTIMATING]" << endl;
      vector<double> yield_estimates;
      vector<double> saturation_estimates;
      double upper_bound, lower_bound = 0.0;
      single_estimates(VERBOSE, counts_hist, orig_max_terms, diagonal, 
		       values_sum, step_size, max_extrapolation, max_val, 
		       val_step, upper_bound, lower_bound,
		       saturation_estimates, yield_estimates);
      if (VERBOSE) 
	cerr << "[WRITING OUTPUT]" << endl;
    
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "#TOTAL_READS\tEXPECTED_DISTINCT" <<  endl;
      double val = 0.0;
      for (size_t i = 0; i < yield_estimates.size(); ++i, val += val_step)
	out << fixed << setprecision(1) << (val + 1.0)*values_sum << '\t' 
	    << yield_estimates[i] << endl;
      
      if (VERBOSE || !size_outfile.empty()) {
	std::ofstream st_of;
	if (!size_outfile.empty()) 
	  st_of.open(size_outfile.c_str());
	ostream size_out(size_outfile.empty() ? cerr.rdbuf() : st_of.rdbuf());
	size_out << "UPPER_BOUND\t" << upper_bound << endl
		 << "LOWER_BOUND\t" << lower_bound << endl;
      }
      
      if (!saturation_outfile.empty()) {
	std::ofstream s_of;
	if (!saturation_outfile.empty())
	  s_of.open(saturation_outfile.c_str());
	ostream sat_out(s_of.rdbuf());
	sat_out << "#TOTAL_REANDS" << "\t"
		<< "SATURATION" << endl;
	val = 0.0;
	for (size_t i = 0; i < saturation_estimates.size(); ++i, val += val_step)
	  out << fixed << setprecision(1) 
	      << (val + 1.0)*values_sum << '\t' 
	      << saturation_estimates[i] << endl;
      }
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
