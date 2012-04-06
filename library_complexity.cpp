/*    library_complexity:
 *
 *    Copyright (C) 2011 University of Southern California and
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <fstream>
#include <numeric>
#include <vector>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::ostream;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;

using smithlab::log_sum_log_vec;

#ifdef HAVE_BAMTOOLS
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;

static SimpleGenomicRegion
BamAlignmentToSimpleGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
				  const BamAlignment &ba) {
  const unordered_map<size_t, string>::const_iterator 
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  
  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = start + ba.Length;
  return SimpleGenomicRegion(chrom, start, end);
}


static void
ReadBAMFormatInput(const string &infile, vector<SimpleGenomicRegion> &read_locations) {
  
  BamReader reader;
  reader.Open(infile);
  
  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  
  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;
  
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
    read_locations.push_back(BamAlignmentToSimpleGenomicRegion(chrom_lookup, bam));
  reader.Close();
}
#endif


static void
get_counts(const vector<SimpleGenomicRegion> &reads, vector<double> &counts) {
  counts.push_back(1);
  for (size_t i = 1; i < reads.size(); ++i)
    if (reads[i] == reads[i - 1]) counts.back()++;
    else counts.push_back(1);
}


static inline double
weight_exponential(const double dist, double decay_factor) {
  return std::pow(0.5, decay_factor*dist);
}

// want #reads in hist_in = reads in hist_out
static void
renormalize_hist(const vector<double> &hist_in,
		 vector<double> &hist_out){
 double out_vals_sum = 0.0;
 for(size_t i = 0; i < hist_out.size(); i++)
   out_vals_sum += i*hist_out[i];
 double in_vals_sum = 0.0;
 for(size_t i = 0; i < hist_in.size(); i++)
   in_vals_sum += i*hist_in[i];
 for(size_t  i =0; i < hist_out.size(); i++)
   hist_out[i] = hist_out[i]*in_vals_sum/out_vals_sum;
}

//additive smoothing
static void
smooth_hist_laplace(const vector<double> &hist_in,
		    const double additive_term,
		    const size_t bandwidth,
		    vector<double> &hist_out){
  size_t hist_out_size = hist_in.size() + bandwidth;
  vector<double> updated_hist(hist_out_size, additive_term);
  updated_hist[0] = 0;
  for(size_t i = 0; i < hist_in.size(); i++)
    updated_hist[i] += hist_in[i];

  renormalize_hist(hist_in, updated_hist);
  hist_out.swap(updated_hist);
}


void
resample_values(const vector<double> &in_values,
		const gsl_rng *rng,
		vector<double> &out_values){
  out_values.clear();
  const double orig_number_captures = accumulate(in_values.begin(), 
						 in_values.end(), 0.0);
  double sampled_number_captures = 0.0;
  while(sampled_number_captures < orig_number_captures){
    double u = gsl_rng_uniform(rng);
    while(u == 1.0)
      u = gsl_rng_uniform(rng);
    const size_t sample_indx  = 
      static_cast<size_t>(floor(in_values.size()*u));
    sampled_number_captures += in_values[sample_indx];
    if(sampled_number_captures > orig_number_captures)
      out_values.push_back(static_cast<double>(sampled_number_captures 
					      - orig_number_captures));
    // ensure that number of captures is constant in resampling
    else
      out_values.push_back(in_values[sample_indx]);
  } 
}

static bool
check_estimates(const vector<double> &estimates) {
  // make sure that the estimate is increasing in the time_step and
  // is below the initial distinct per step_size
  if(!finite(accumulate(estimates.begin(), estimates.end(), 0.0)) || 
     estimates.size() == 0)
    return false;

    for (size_t i = 1; i < estimates.size(); ++i)
    if (estimates[i] < estimates[i - 1] ||
	(i >= 2 && (estimates[i] - estimates[i - 1] >
		    estimates[i - 1] - estimates[i - 2])) ||
	estimates[i] < 0.0)
      return false;
  

  return true;
}

void
laplace_bootstrap_smoothed_hist(const bool VERBOSE, const vector<double> &orig_values, 
				const double smoothing_val,
				const size_t bootstraps, const size_t orig_max_terms,
				const int diagonal, const double step_size,
				const double max_extrapolation, const double max_val,
				const double val_step, const size_t bandwidth,
				vector<double> &lower_bound_size,
				vector<double> &upper_bound_size,
				vector< vector<double> > &lower_estimates){
  // clear returning vectors
  upper_bound_size.clear();
  lower_bound_size.clear();
  lower_estimates.clear();

  //setup rng
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  const int seed = time(0) + getpid();
  srand(seed);
  gsl_rng_set(rng, rand()); 
  

  if(VERBOSE) cerr << "lower" << endl;
  int iter = 0;
  while(lower_estimates.size() < bootstraps){
    if(VERBOSE) cerr << iter << "\t" << lower_estimates.size() << endl;
    iter++;

    vector<double> lower_boot_estimates;

    bool ACCEPT_ESTIMATES = false;

    vector<double> boot_values;

    //  resample_histogram(rng, smooth_hist, boot_hist);

    resample_values(orig_values, rng, boot_values);

    const size_t max_observed_count = 
      *std::max_element(boot_values.begin(), boot_values.end());

    // BUILD THE HISTOGRAM
    vector<double> boot_hist(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < boot_values.size(); ++i)
      ++boot_hist[boot_values[i]];


    //resize boot_hist to remove excess zeros
    while(boot_hist.back() == 0)
      boot_hist.pop_back();

    vector<double> smooth_boot_hist;
    smooth_hist_laplace(boot_hist, smoothing_val, bandwidth, smooth_boot_hist);

    if(VERBOSE){
      cerr << "boot_hist \t" << boot_hist.size() << endl;
      for(size_t i = 0; i < boot_hist.size(); i++)
	cerr << boot_hist[i] << "\t";
      cerr << endl;

      cerr << "smoothed boot_hist \t" << smooth_boot_hist.size() << endl;
      for(size_t i = 0; i < smooth_boot_hist.size(); i++)
	cerr << smooth_boot_hist[i] << "\t";
      cerr << endl;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < smooth_boot_hist.size()  && 
	   smooth_boot_hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);  

   //refit curve for lower bound
    max_terms = max_terms - (max_terms % 2 == 1);

    //refit curve for lower bound
    const ContinuedFractionApproximation lower_cfa(diagonal, max_terms, 
						     step_size, max_extrapolation);
    const ContinuedFraction lower_cf(lower_cfa.optimal_continued_fraction(smooth_boot_hist));

    lower_boot_estimates.clear();
    if(lower_cf.is_valid())
      lower_cf.extrapolate_distinct(smooth_boot_hist, max_val, 
				    val_step, lower_boot_estimates);
      
      //sanity check
    ACCEPT_ESTIMATES = check_estimates(lower_boot_estimates);
    if(VERBOSE) cerr << ACCEPT_ESTIMATES << endl;
    if(ACCEPT_ESTIMATES){
      const double distinct  = accumulate(boot_hist.begin(), 
					  boot_hist.end(), 0.0);
      const double upper_bound =
	upperbound_librarysize(boot_hist, lower_cf.return_degree()) + distinct;
      const double lower_bound = 
	lower_cfa.lowerbound_librarysize(VERBOSE, smooth_boot_hist, upper_bound);
      if(finite(lower_bound) && lower_bound > 0 && finite(upper_bound)){
	lower_estimates.push_back(lower_boot_estimates);
	upper_bound_size.push_back(upper_bound);
	lower_bound_size.push_back(lower_bound);
	if(VERBOSE)
	  cerr << "upper_bound \t" << upper_bound_size.back() << endl;
      }
    }
  }
}

static inline double
compute_var(const vector<double> &estimates){
  const double sample_size = static_cast<double>(estimates.size());
  double mean = 
    accumulate(estimates.begin(), estimates.end(), 0.0)/sample_size;
  double variance = 0.0;
  for(size_t i = 0; i < estimates.size(); i++)
    variance += (estimates[i] - mean)*(estimates[i] - mean)/sample_size;
  return variance; 
}


static inline double
alpha_log_confint_multiplier(const double estimate,
			     const double initial_distinct,
			     const double variance,
			     const double alpha){
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  return exp(inv_norm_alpha*sqrt(log(1.0 + 
				     variance/pow(estimate - 
						  initial_distinct, 2))));
}

static void
return_median_and_alphaCI(const vector< vector<double> > &estimates,
			  const double alpha, const double initial_distinct,
			  vector<double> &median_estimates,
			  vector<double> &lower_alphaCI,
			  vector<double> &upper_alphaCI){

  for(size_t i = 0; i < estimates[0].size(); i++){
    // estimates is in wrong order, work locally on const val
    vector<double> estimates_row(estimates.size(), 0.0);
    for(size_t k = 0; k < estimates_row.size(); ++k)
      estimates_row[k] = estimates[k][i];
      
    //   const double mean = compute_mean(estimates_row);
    //   mean_estimates.push_back(mean);
    //sort to get confidence interval
    sort(estimates_row.begin(), estimates_row.end());
    median_estimates.push_back(estimates_row[estimates.size()/2 - 1]); //median
    const double variance = compute_var(estimates_row);
    const double confint_multiplier = 
      alpha_log_confint_multiplier(median_estimates.back(), initial_distinct,
				   variance, alpha);
    upper_alphaCI.push_back(initial_distinct + 
			    (median_estimates.back() - initial_distinct)*confint_multiplier);
    lower_alphaCI.push_back(initial_distinct + 
			    (median_estimates.back() - initial_distinct)/confint_multiplier);		    
  }

}






int
main(const int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    string library_size_outfile;

    size_t orig_max_terms = 100;
    double max_extrapolation = 1e10;
    double step_size = 1e6;
    size_t smoothing_bandwidth = 4;
    double smoothing_val = 1e-3;
    size_t bootstraps = 100;
    int diagonal = -1;
    double alpha = 0.05;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool SMOOTH_HISTOGRAM = false;  
    
#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(argv[0], "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("LIBRARY_SIZE", 'L', "library size output file", 
		      false , library_size_outfile);
    opt_parse.add_opt("extrapolation_length",'e',
		      "maximum extrapolation length, default 1e10", 
                      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size between extrapolations, default 1e6", 
                      false, step_size);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps, default 100",
		      false, bootstraps);
    opt_parse.add_opt("alpha", 'a', "alpha for confidence intervals, default 0.05",
		      false, alpha);
    opt_parse.add_opt("terms",'t',"maximum number of terms", false, orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
#ifdef HAVE_BAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format", 
		      false , BAM_FORMAT_INPUT);
#endif

    opt_parse.add_opt("smooth",'\0',"smooth histogram (default: no smoothing)",
		      false, SMOOTH_HISTOGRAM);


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
    /**********************************************************************/
    
    // READ IN THE DATA
    vector<SimpleGenomicRegion> read_locations;
#ifdef HAVE_BAMTOOLS
    if (BAM_FORMAT_INPUT)
      ReadBAMFormatInput(input_file_name, read_locations);
    else 
#endif
      ReadBEDFile(input_file_name, read_locations);
    if (!check_sorted(read_locations))
      throw SMITHLABException("read_locations not sorted");

    // OBTAIN THE COUNTS FOR DISTINCT READS
    vector<double> values;
    get_counts(read_locations, values);
    const double initial_distinct = static_cast<double>(values.size());
    
    // JUST A SANITY CHECK
    const size_t n_reads = read_locations.size();
    const size_t values_sum = accumulate(values.begin(), values.end(), 0ul);
    assert(values_sum == n_reads);
    
    const double max_val = max_extrapolation/static_cast<double>(values_sum);
    const double val_step = step_size/static_cast<double>(values_sum);
    
    const size_t max_observed_count = 
      *std::max_element(values.begin(), values.end());

    // BUILD THE HISTOGRAM
    vector<double> counts_histogram(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_histogram[values[i]];
    
    const size_t distinct_counts = 
      std::count_if(counts_histogram.begin(), counts_histogram.end(),
		    bind2nd(std::greater<size_t>(), 0));
    
    if (VERBOSE)
      cerr << "TOTAL READS     = " << read_locations.size() << endl
	   << "DISTINCT COUNTS = " << distinct_counts << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_histogram[1] << endl
	   << "MAX TERMS       = " << orig_max_terms << endl;


    if(VERBOSE){
      cerr << "orig_hist \t" << counts_histogram.size() << endl;
      for(size_t i = 0; i < counts_histogram.size(); i++)
	cerr << counts_histogram[i] << "\t";
      cerr << endl;
    }
    

    if(!SMOOTH_HISTOGRAM){
      smoothing_val = 0.0;
      smoothing_bandwidth = 0;
    }

  
    if(VERBOSE) cerr << "laplace resampling" << endl;
    vector<vector <double> > boot_estimates;
    vector<double> lower_librarysize;
    vector<double> upper_librarysize;
    laplace_bootstrap_smoothed_hist(VERBOSE, values, smoothing_val, bootstraps, orig_max_terms,
				    diagonal, step_size, max_extrapolation, max_val,
				    val_step, smoothing_bandwidth, 
				    lower_librarysize, upper_librarysize,
				    boot_estimates);

    if(VERBOSE) cerr << "compute confidence intervals" << endl;

    vector<double> median_estimates;
    vector<double> upper_alphaCI;
    vector<double> lower_alphaCI;
    return_median_and_alphaCI(boot_estimates, alpha, initial_distinct,
			      median_estimates, lower_alphaCI, upper_alphaCI);

    if(VERBOSE) cerr << "outputing" << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    double val = 0.0;
    out << "reads" << '\t' << "median_estimates" << '\t' 
	<< "lowerCI" << '\t' << "upperCI" << "\t" 
	<< endl;
    for (size_t i = 0; i < median_estimates.size(); ++i, val += val_step)
      out << std::fixed << std::setprecision(1) 
	  << (val + 1.0)*values_sum << '\t' << median_estimates[i] << '\t'
	  << lower_alphaCI[i] << '\t' << upper_alphaCI[i] << endl;

    if (VERBOSE || !library_size_outfile.empty()) {
      std::ofstream stats_of;
      if (!library_size_outfile.empty()) stats_of.open(library_size_outfile.c_str());
      ostream stats_out(library_size_outfile.empty() ? cerr.rdbuf() : stats_of.rdbuf());
      const double cf_upper_size_var = compute_var(upper_librarysize);
      const double cf_upper_size_mean = 
	accumulate(upper_librarysize.begin(), upper_librarysize.end(), 0.0)/upper_librarysize.size();
      const double cf_upper_size_alpha_multiplier = 
	alpha_log_confint_multiplier(cf_upper_size_mean, initial_distinct, 
				     cf_upper_size_var, alpha);
      stats_out << "LIBRARY_SIZE_UPPERBOUND_MEAN\t" << cf_upper_size_mean << endl;
      stats_out << "LIBRARY_SIZE_UPPERBOUND_LOWER" << 100*(1-alpha) << "%CI\t"
		<< initial_distinct + (cf_upper_size_mean - initial_distinct)/cf_upper_size_alpha_multiplier
		<< endl;
      stats_out << "LIBRARY_SIZE_UPPERBOUND_UPPER" << 100*(1-alpha) << "%CI\t"
		<< initial_distinct + (cf_upper_size_mean - initial_distinct)*cf_upper_size_alpha_multiplier
		<< endl;
      const double cf_lower_size_var = compute_var(lower_librarysize);
      const double cf_lower_size_mean = 
	accumulate(lower_librarysize.begin(), lower_librarysize.end(), 0.0)/lower_librarysize.size();
      const double cf_lower_size_alpha_multiplier = 
	alpha_log_confint_multiplier(cf_lower_size_mean, initial_distinct, 
				     cf_lower_size_var, alpha);
      stats_out << "LIBRARY_SIZE_LOWERBOUND_MEAN\t" << cf_lower_size_mean << endl;
      stats_out << "LIBRARY_SIZE_LOWERBOUND_LOWER" << 100*(1-alpha) << "%CI\t"
		<< initial_distinct + (cf_lower_size_mean - initial_distinct)/cf_lower_size_alpha_multiplier
		<< endl;
      stats_out << "LIBRARY_SIZE_LOWERBOUND_UPPER" << 100*(1-alpha) << "%CI\t"
		<< initial_distinct + (cf_lower_size_mean - initial_distinct)*cf_lower_size_alpha_multiplier
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
