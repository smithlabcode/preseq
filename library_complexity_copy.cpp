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
#include "ZTP.hpp"

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
get_counts(const vector<SimpleGenomicRegion> &reads, vector<size_t> &counts) {
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

// exponential smoothing
static void
smooth_histogram(const size_t bandwidth,
                 const double decay_factor, vector<double> &hist) {
  vector<double> updated_hist(hist);
  for (size_t i = 0; i < hist.size(); ++i) {
    double total_prob = 0.0, total_weight = 0.0;
    for (size_t j = ((i >= bandwidth/2) ? i - bandwidth/2 : 0);
         j < std::min(i + bandwidth/2, hist.size()); ++j) {
      const double dist = std::abs(int(i) - int(j));
      const double w = weight_exponential(dist, decay_factor);
      total_prob += w*hist[j];
      total_weight += w;
    }
    updated_hist[i] = total_prob/total_weight;
  }

  //renormalize hist so that number of captures is constant
  updated_hist[0] = 0;
  renormalize_hist(hist, updated_hist);

  updated_hist.swap(hist);
}

//additive smoothing
static void
smooth_hist_laplace(const vector<double> &hist_in,
		    const double additive_term,
		    const size_t bandwidth,
		    vector<double> &hist_out){
  hist_out.clear();
  size_t hist_out_size = hist_in.size() + bandwidth;
  hist_out.resize(hist_out_size, additive_term);
  hist_out[0] = 0;
  for(size_t i = 0; i < hist_in.size(); i++)
    hist_out[i] += hist_in[i];
  renormalize_hist(hist_in, hist_out);
}


// to bootstrap histogram
static void
resample_histogram(const gsl_rng *rng,
		   const vector<double> &hist_in,
		   vector<double> &hist_out) {

 const double n_samples = std::accumulate(hist_in.begin(), hist_in.end(), 0.0);

 vector<double> cumulants(1, 0.0);
 for (size_t i = 0; i < hist_in.size(); ++i)
   cumulants.push_back(cumulants.back() + hist_in[i]);

 hist_out.resize(hist_in.size(), 0.0);
 for (size_t i = 0; i < n_samples; ++i) {
   const size_t idx = std::lower_bound(cumulants.begin(), cumulants.end(), 
				       n_samples*gsl_rng_uniform(rng)) - cumulants.begin();
   hist_out[idx-1]++;
 }

 //renormalize hist_out so that number of captures is same
 renormalize_hist(hist_in, hist_out);
}


// smooth histogram a local ZTP, bin width ~ log(val)
static void
smooth_histogram(const vector<double> &hist_in,
		 const double bin_radius_multiplier, const size_t start_indx,
		 vector<double> &hist_out){
  assert(hist_in.size() > 0);
  hist_out.clear();
  //new hist should be larger in size since we are smoothing max vals in both directions
  vector<double> weighted_hist(hist_in.size() + 
			       static_cast<size_t>(ceil(bin_radius_multiplier*log(hist_in.size()))), 0.0);


  //count times so we average the probabilities in each cell seperately
  vector<double> count_times_visited(weighted_hist);
  size_t lower_bound = start_indx;
  

  for(size_t i = start_indx; i < hist_in.size(); i++){
    const size_t bin_radius = static_cast<size_t>(floor(bin_radius_multiplier*log(i)));
    if(bin_radius > 0){
      const size_t iter_lower_bound = std::max(static_cast<size_t>(ceil(i - bin_radius)), lower_bound);
      const size_t iter_upper_bound = std::min(i + bin_radius + 2, weighted_hist.size());

      vector<double> local_histogram(iter_upper_bound+1, 0.0);
      for(size_t j = iter_lower_bound; j < std::min(iter_upper_bound,
						     hist_in.size()); ++j){
	local_histogram[j] += hist_in[j];
      }

    // average weight in window
      const double local_weight = accumulate(local_histogram.begin(), local_histogram.end(), 0.0);

      if(local_weight > 0.0){ // smooth only if there is weight in window
	   //fit ZTP in window
		ZeroTruncatedPoisson local_distro(1.0);
	local_distro.estimateParams(local_histogram, false);

      // only add weight if in window proportional to the likelihood
	vector<double> local_probs;
	for(size_t k = iter_lower_bound; k <= iter_upper_bound; k++)
	  local_probs.push_back(exp(local_distro.logLikelihood(k)));
       
	for(size_t j = iter_lower_bound; j < iter_upper_bound; ++j){
	  	  weighted_hist[j] += 
	  	    local_weight*local_probs[j - iter_lower_bound];
       	  count_times_visited[j]++;
	}
      }
    }
    else{
      weighted_hist[i] += hist_in[i];
      count_times_visited[i]++;
      lower_bound++;
    }
  
  }
  //renormalize so total weight is constant and we are averaging the weight
  //to each pos
  double total_weighted_captures = 0.0;
  for(size_t i = start_indx; i < weighted_hist.size(); i++){
    if(count_times_visited[i] > 0)
      weighted_hist[i] = weighted_hist[i]/count_times_visited[i];
    total_weighted_captures += i*weighted_hist[i];
  }

  double total_unweighted_captures = 0.0;
  for(size_t i = start_indx; i < hist_in.size(); i++)
    total_unweighted_captures += i*hist_in[i];

  for(size_t  i = start_indx; i < weighted_hist.size(); i++)
    weighted_hist[i] = weighted_hist[i]*total_unweighted_captures/total_weighted_captures;

  for(size_t  i = 0; i < start_indx + 1; i++)
    weighted_hist[i] = hist_in[i];

  hist_out.swap(weighted_hist);
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
laplace_bootstrap_smoothed_hist(const bool VERBOSE, const vector<double> &smooth_hist, const double additive_val,
			    const size_t bootstraps, const size_t orig_max_terms,
			    const size_t diagonal, const double step_size,
			    const double max_extrapolation, const double max_val,
				const double val_step, const size_t bandwidth,
			    vector< vector< double> > &lower_estimates){

  //setup rng
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  const int seed = time(0) + getpid();
  srand(seed);
  gsl_rng_set(rng, rand()); 

  for(size_t i = 0; i < bootstraps; i++){
    if(VERBOSE) cerr << i << endl;
    vector<double> lower_boot_estimates;

    bool ACCEPT_ESTIMATES = false;
    while(!ACCEPT_ESTIMATES){

      vector<double> boot_hist(smooth_hist);
      resample_histogram(rng, smooth_hist, boot_hist);

      vector<double> smooth_boot_hist;
      smooth_hist_laplace(boot_hist, additive_val, bandwidth, smooth_boot_hist);

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
	lower_cf.extrapolate_distinct(VERBOSE, smooth_boot_hist, max_val, 
				      val_step, lower_boot_estimates);

      //sanity check
      ACCEPT_ESTIMATES = check_estimates(lower_boot_estimates);
      if(VERBOSE) cerr << ACCEPT_ESTIMATES << endl;
    }
    lower_estimates.push_back(lower_boot_estimates);
  }
}


void
ztp_bootstrap_smoothed_hist(const bool VERBOSE, const vector<double> &smooth_hist,
			    const double bin_radius, const size_t start_indx,  
			    const size_t bootstraps, const size_t orig_max_terms,
			    const size_t diagonal, const double step_size,
			    const double max_extrapolation, const double max_val,
			    const double val_step,
			    vector<vector< double> > &lower_estimates){
  //setup rng
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  srand(time(0) + getpid());
  gsl_rng_set(rng, rand()); 
  for(size_t i = 0; i < bootstraps; i++){

    vector<double> lower_boot_estimates;
    bool ACCEPT_ESTIMATES = false;
    while(!ACCEPT_ESTIMATES){

      vector<double> boot_hist(smooth_hist);
      resample_histogram(rng, smooth_hist, boot_hist);

      vector<double> smooth_boot_hist(boot_hist);
      smooth_histogram(boot_hist, bin_radius, start_indx, smooth_boot_hist);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
      size_t counts_before_first_zero = 1;
      while (counts_before_first_zero < smooth_boot_hist.size() && 
	     smooth_boot_hist[counts_before_first_zero] > 0)
	++counts_before_first_zero;
      size_t max_terms = std::min(orig_max_terms, counts_before_first_zero);  

   //refit curve for lower bound
      max_terms = max_terms - (max_terms % 2 == 1);

    //refit curve for lower bound
      const ContinuedFractionApproximation lower_cfa(diagonal, max_terms, 
						     step_size, max_extrapolation);
      const ContinuedFraction lower_cf(lower_cfa.optimal_continued_fraction(smooth_boot_hist));

      lower_boot_estimates.clear();
      if(lower_cf.is_valid())
	lower_cf.extrapolate_distinct(VERBOSE, smooth_boot_hist, max_val, 
				      val_step, lower_boot_estimates);

 
      //sanity check
      ACCEPT_ESTIMATES = check_estimates(lower_boot_estimates);
    }
    lower_estimates.push_back(lower_boot_estimates);
  }
}

void
andrew_bootstrap_smoothed_hist(const bool VERBOSE, const vector<double> &smooth_hist,
			       const size_t bootstraps, const size_t orig_max_terms, 
			       const size_t bandwidth, const double decay_factor,
			       const size_t diagonal, const double step_size,
			       const double max_extrapolation, const double max_val,
			       const double val_step,
			       vector< vector<double> > &lower_estimates){
  //setup rng
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  srand(time(0) + getpid());
  gsl_rng_set(rng, rand()); 
  for(size_t i = 0; i < bootstraps; i++){
    vector<double> lower_boot_estimates;
    bool ACCEPT_ESTIMATES = false;
    while(!ACCEPT_ESTIMATES){

      vector<double> boot_hist(smooth_hist);
      resample_histogram(rng, smooth_hist, boot_hist);

      smooth_histogram(bandwidth, decay_factor, boot_hist);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
      size_t counts_before_first_zero = 1;
      while (counts_before_first_zero < boot_hist.size() && 
	     boot_hist[counts_before_first_zero] > 0)
	++counts_before_first_zero;
      size_t max_terms = std::min(orig_max_terms, counts_before_first_zero);  

      //ensure we are using underestimate
      max_terms = max_terms - (max_terms % 2 == 1);

    //refit curve for lower bound
      const ContinuedFractionApproximation lower_cfa(diagonal, max_terms, 
						     step_size, max_extrapolation);
      const ContinuedFraction lower_cf(lower_cfa.optimal_continued_fraction(boot_hist));

      lower_boot_estimates.clear();
      if(lower_cf.is_valid())
	lower_cf.extrapolate_distinct(VERBOSE, boot_hist, max_val, 
				      val_step, lower_boot_estimates);

      //sanity check
      ACCEPT_ESTIMATES = 
	check_estimates(lower_boot_estimates);
    }
    lower_estimates.push_back(lower_boot_estimates);
  }
}

void
bootstrap_smoothed_hist(const bool VERBOSE, const vector<double> &smooth_hist,
			const size_t bootstraps, const size_t orig_max_terms,
			const size_t diagonal, const double step_size,
			const double max_extrapolation, const double max_val,
			const double val_step,
			vector< vector<double> > &lower_estimates){
  //setup rng
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  srand(time(0) + getpid());
  gsl_rng_set(rng, rand()); 
  for(size_t i = 0; i < bootstraps; i++){
    vector<double> lower_boot_estimates;
    bool ACCEPT_ESTIMATES = false;
    while(!ACCEPT_ESTIMATES){
      vector<double> boot_hist(smooth_hist);
      resample_histogram(rng, smooth_hist, boot_hist);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
      size_t counts_before_first_zero = 1;
      while (counts_before_first_zero < boot_hist.size() && 
	     boot_hist[counts_before_first_zero] > 0)
	++counts_before_first_zero;
      size_t max_terms = std::min(orig_max_terms, counts_before_first_zero);  

    //refit curve for lower bound
      max_terms = max_terms - (max_terms % 2 == 1);

     //refit curve for lower bound
      const ContinuedFractionApproximation lower_cfa(diagonal, max_terms, 
						     step_size, max_extrapolation);
      const ContinuedFraction lower_cf(lower_cfa.optimal_continued_fraction(boot_hist));

      lower_boot_estimates.clear();
      if(lower_cf.is_valid())
	lower_cf.extrapolate_distinct(VERBOSE, boot_hist, max_val, 
				      val_step, lower_boot_estimates);

      //sanity check
      ACCEPT_ESTIMATES = 
	check_estimates(lower_boot_estimates);
    }
    lower_estimates.push_back(lower_boot_estimates);    
  }
}

static inline double
compute_mean(const vector<double> &estimates){
  return accumulate(estimates.begin(), estimates.end(), 0.0)/estimates.size();
}

void
return_mean_and_alphaCI(const vector< vector<double> > &estimates,
			const double alpha, 
			vector<double> &mean_estimates,
			vector<double> &lower_CI, 
			vector<double> &upper_CI){
  const size_t lower_alpha_percentile = 
    static_cast<size_t>(floor(alpha*estimates.size()/2));

  const size_t upper_alpha_percentile = estimates.size() - lower_alpha_percentile;

  for(size_t i = 0; i < estimates[0].size(); i++){
    // estimates is in wrong order, work locally on const val
    vector<double> estimates_row(estimates.size(), 0.0);
    for(size_t k = 0; k < estimates_row.size(); ++k)
      estimates_row[k] = estimates[k][i];
      
    const double mean = compute_mean(estimates_row);
    mean_estimates.push_back(mean);
    //sort to get confidence interval
    sort(estimates_row.begin(), estimates_row.end());
    lower_CI.push_back(estimates_row[lower_alpha_percentile]);
    upper_CI.push_back(estimates_row[upper_alpha_percentile]);
  }

}






int
main(const int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    string stats_outfile;

    size_t orig_max_terms = 100;
    double max_extrapolation = 1e10;
    double step_size = 1e6;
    size_t smoothing_bandwidth = 4;
    double smoothing_decay_factor = 15.0;
    double smooth_val = 1e-4;
    size_t bootstraps = 1000;
    int diagonal = -1;
    double bin_radius = 1.0;
    size_t start_indx = 1;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool SMOOTH_HISTOGRAM = true; //false; 
    // bool LIBRARY_SIZE = false;
    
#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(argv[0], "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("stats", 'S', "stats output file", 
		      false , stats_outfile);
    opt_parse.add_opt("extrapolation_length",'e',"maximum extrapolation length", 
                      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size between extrapolations", 
                      false, step_size);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps",
		      false, bootstraps);
    opt_parse.add_opt("smooth_val",'m',"value to smooth by in additive",
		      false, smooth_val);
    opt_parse.add_opt("smoothing_bandwidth",'w'," ",
		      false, smoothing_bandwidth);
    opt_parse.add_opt("decay_factor",'c',"smoothing_ decay_factor",
		      false, smoothing_decay_factor);
    opt_parse.add_opt("bin_radius",'r',"bin size radius multiplier",
		      false, bin_radius);
    opt_parse.add_opt("start_indx",'x',"start index of smoothing",
		      false, start_indx);
    opt_parse.add_opt("diagonal",'d',"diagonal for cont frac",
		      false, diagonal);
    opt_parse.add_opt("terms",'t',"maximum number of terms", false, orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
#ifdef HAVE_BAMTOOLS
    opt_parse.add_opt("bam", 'b', "input is in BAM format", 
		      false , BAM_FORMAT_INPUT);
#endif
    //     opt_parse.add_opt("tol", '\0', "general numerical tolerance",
    // 		      false, tolerance);
    //     opt_parse.add_opt("delta", '\0', "derivative step size",
    //                       false, deriv_delta);
    //     opt_parse.add_opt("smooth",'\0',"smooth histogram (default: no smoothing)",
    //                       false, SMOOTH_HISTOGRAM);
    //     opt_parse.add_opt("bandwidth", '\0', "smoothing bandwidth",
    // 		      false, smoothing_bandwidth);
    //     opt_parse.add_opt("decay", '\0', "smoothing decay factor",
    // 		      false, smoothing_decay_factor);
    //     opt_parse.add_opt("diag", 'd', "diagonal to use",
    // 		      false, diagonal);
    //     opt_parse.add_opt("library_size", '\0', "estimate library size "
    // 		      "(default: estimate distinct)", false, LIBRARY_SIZE);
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
    vector<size_t> values;
    get_counts(read_locations, values);
    
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

    vector<double> andrew_counts_hist(counts_histogram);
    vector<double> ztp_counts_hist(counts_histogram);
    vector<double> laplace_counts_hist(counts_histogram);
    if (SMOOTH_HISTOGRAM){ 
      smooth_histogram(smoothing_bandwidth, 
		       smoothing_decay_factor, andrew_counts_hist);
      smooth_histogram(counts_histogram, bin_radius, start_indx, ztp_counts_hist);
      smooth_hist_laplace(counts_histogram, smooth_val, smoothing_bandwidth, laplace_counts_hist);
    }

    if(VERBOSE){
      cerr << "orig_hist" << endl;
      for(size_t i = 0; i < counts_histogram.size(); i++)
	cerr << counts_histogram[i] << "\t";
      cerr << endl;

      cerr << "smoothed_hist_andrew" << endl;
      for(size_t i = 0; i < andrew_counts_hist.size(); i++)
	cerr << andrew_counts_hist[i] << "\t";
      cerr << endl;

      cerr << "smoothed_hist_ZTP" << endl;
      for(size_t i = 0; i < ztp_counts_hist.size(); i++)
	cerr << ztp_counts_hist[i] << "\t";
      cerr << endl;

      cerr << "laplace_smoothed_hist" << endl;
      for(size_t i = 0; i < laplace_counts_hist.size(); i++)
	cerr << laplace_counts_hist[i] << "\t";
      cerr << endl;
    }

    cerr << "laplace" << endl;
    vector<vector <double> > lower_laplace_boot_estimates;
    laplace_bootstrap_smoothed_hist(VERBOSE, laplace_counts_hist, smooth_val, bootstraps, orig_max_terms,
				    diagonal, step_size, max_extrapolation, max_val,
				    val_step,smoothing_bandwidth, 
				    lower_laplace_boot_estimates);


    /*
    cerr << "andrew" << endl;
    vector< vector<double> > lower_andrew_boot_estimates;
    vector< vector<double> > upper_andrew_boot_estimates;
    bootstrap_smoothed_hist(andrew_counts_hist, bootstraps, 
			    orig_max_terms, diagonal, step_size,
			    max_extrapolation, max_val, val_step,
			    lower_andrew_boot_estimates,
			    upper_andrew_boot_estimates);
    cerr << "ztp " << ztp_counts_hist.size() << endl;
    vector< vector<double> > lower_ztp_boot_estimates;
    vector< vector<double> > upper_ztp_boot_estimates;
    bootstrap_smoothed_hist(ztp_counts_hist, bootstraps, 
			    orig_max_terms, diagonal, step_size,
			    max_extrapolation, max_val, val_step,
			    lower_ztp_boot_estimates,
			    upper_ztp_boot_estimates);
    cerr << "andrew_smooth" << endl;
    vector< vector<double> > lower_andrew_smooth_boot_estimates;
    vector< vector<double> > upper_andrew_smooth_boot_estimates;
    andrew_bootstrap_smoothed_hist(andrew_counts_hist, bootstraps,
				   orig_max_terms, smoothing_bandwidth,
				   smoothing_decay_factor, diagonal, step_size,
			    max_extrapolation, max_val, val_step, lower_andrew_smooth_boot_estimates,
				   upper_andrew_smooth_boot_estimates);
    cerr << "ztp_smooth" << endl;
    vector< vector<double> > lower_ztp_smooth_boot_estimates;
    vector< vector<double> > upper_ztp_smooth_boot_estimates;
    ztp_bootstrap_smoothed_hist(ztp_counts_hist, bootstraps, orig_max_terms,
				bin_radius, start_indx,  diagonal, step_size,
			    max_extrapolation, max_val, val_step,
				lower_ztp_smooth_boot_estimates,
				upper_ztp_smooth_boot_estimates);
    

    vector<double> lower_andrew_boot_mean;
    vector<double> lower_andrew_boot_lowerCI;
    vector<double> lower_andrew_boot_upperCI;
    return_mean_and_alphaCI(lower_andrew_boot_estimates, 0.05, lower_andrew_boot_mean,
			      lower_andrew_boot_lowerCI, lower_andrew_boot_upperCI);




    vector<double> upper_andrew_boot_mean;
    vector<double> upper_andrew_boot_lowerCI;
    vector<double> upper_andrew_boot_upperCI;
    return_mean_and_alphaCI(upper_andrew_boot_estimates, 0.05, upper_andrew_boot_mean,
			      upper_andrew_boot_lowerCI, upper_andrew_boot_upperCI);



    vector<double> lower_ztp_boot_mean;
    vector<double> lower_ztp_boot_lowerCI;
    vector<double> lower_ztp_boot_upperCI;
    return_mean_and_alphaCI(lower_ztp_boot_estimates, 0.05, lower_ztp_boot_mean,
			      lower_ztp_boot_lowerCI, lower_ztp_boot_upperCI);


    vector<double> upper_ztp_boot_mean;
    vector<double> upper_ztp_boot_lowerCI;
    vector<double> upper_ztp_boot_upperCI;
    return_mean_and_alphaCI(upper_ztp_boot_estimates, 0.05, upper_ztp_boot_mean,
			      upper_ztp_boot_lowerCI, upper_ztp_boot_upperCI);

 
    
   vector<double> lower_andrew_smooth_boot_mean;
    vector<double> lower_andrew_smooth_boot_lowerCI;
    vector<double> lower_andrew_smooth_boot_upperCI;
    return_mean_and_alphaCI(lower_andrew_smooth_boot_estimates, 0.05, lower_andrew_smooth_boot_mean,
			      lower_andrew_smooth_boot_lowerCI, lower_andrew_smooth_boot_upperCI);


    vector<double> upper_andrew_smooth_boot_mean;
    vector<double> upper_andrew_smooth_boot_lowerCI;
    vector<double> upper_andrew_smooth_boot_upperCI;
    return_mean_and_alphaCI(upper_andrew_smooth_boot_estimates, 0.05, upper_andrew_smooth_boot_mean,
			      upper_andrew_smooth_boot_lowerCI, upper_andrew_smooth_boot_upperCI);


    vector<double> lower_ztp_smooth_boot_mean;
    vector<double> lower_ztp_smooth_boot_lowerCI;
    vector<double> lower_ztp_smooth_boot_upperCI;
    return_mean_and_alphaCI(lower_ztp_smooth_boot_estimates, 0.05, lower_ztp_smooth_boot_mean,
			      lower_ztp_smooth_boot_lowerCI, lower_ztp_smooth_boot_upperCI);


    vector<double> upper_ztp_smooth_boot_mean;
    vector<double> upper_ztp_smooth_boot_lowerCI;
    vector<double> upper_ztp_smooth_boot_upperCI;
    return_mean_and_alphaCI(upper_ztp_smooth_boot_estimates, 0.05, upper_ztp_smooth_boot_mean,
			      upper_ztp_smooth_boot_lowerCI, upper_ztp_smooth_boot_upperCI);
    */

    cerr << "compute mean" << endl;

    vector<double> lower_laplace_smooth_boot_mean;
    vector<double> lower_laplace_smooth_boot_lowerCI;
    vector<double> lower_laplace_smooth_boot_upperCI;
    return_mean_and_alphaCI(lower_laplace_boot_estimates, 0.05, lower_laplace_smooth_boot_mean,
			    lower_laplace_smooth_boot_lowerCI, lower_laplace_smooth_boot_upperCI);

    cerr << "outputing" << endl;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    double val = 0.0;
    out << "reads" << '\t' << "lower_laplace_mean" << '\t' 
      << "lower_laplace_lowerCI" << '\t' << "lower_laplace_upper_CI" 
      /* << '\t'
	<< "lower_andrew_mean" << '\t' << "lower_andrew_lowerCI"
	  << '\t' << "lower_andrew_upperCI" << '\t' << "lower_ztp_mean" << '\t' << "lower_ztp_lowerCI"
	  << '\t' << "lower_ztp_upperCI" << '\t' << "lower_andrew_smooth_mean" << '\t' << "lower_andrew_smooth_lowerCI"
	  << '\t' << "andrew_smooth_upperCI" << '\t' << "lower_ztp_smooth_mean"
	  << '\t' << "lower_ztp_smooth_lowerCI"
	  << '\t' << "lower_ztp_smooth_upperCI" */
	<< endl;
    for (size_t i = 0; i < lower_laplace_smooth_boot_mean.size(); ++i, val += val_step)
      out << std::fixed << std::setprecision(1) 
	  << (val + 1.0)*values_sum << '\t' << lower_laplace_smooth_boot_mean[i] << '\t'
	  << lower_laplace_smooth_boot_lowerCI[i] << '\t' << lower_laplace_smooth_boot_upperCI[i] << '\t' 
	/*	lower_andrew_boot_mean[i] << '\t' << lower_andrew_boot_lowerCI[i]
	  << '\t' << lower_andrew_boot_upperCI[i] << '\t' << lower_ztp_boot_mean[i] << '\t' << lower_ztp_boot_lowerCI[i]
	  << '\t' << lower_ztp_boot_upperCI[i] << '\t' << lower_andrew_smooth_boot_mean[i] << '\t' << lower_andrew_smooth_boot_lowerCI[i]
	  << '\t' << lower_andrew_smooth_boot_upperCI[i] << '\t' << lower_ztp_smooth_boot_mean[i] 
	  << '\t' << lower_ztp_smooth_boot_lowerCI[i]
	  << '\t' << lower_ztp_smooth_boot_upperCI[i] */
	  << endl;
    
    /* Finally output the bounds on the library size if either verbose
       output is requested or this info was specifically requested.
     */
    /*   if (VERBOSE || !stats_outfile.empty()) {
      std::ofstream stats_of;
      if (!stats_outfile.empty()) stats_of.open(stats_outfile.c_str());
      ostream stats_out(stats_outfile.empty() ? cerr.rdbuf() : stats_of.rdbuf());
      
      const double upper_bound = distinct_counts +
	upperbound_librarysize(counts_histogram, max_terms);
      stats_out << "CF_UPPER=" << upper_bound << endl;
      stats_out << "CF_LOWER="
		<< cfa.lowerbound_librarysize(counts_histogram, upper_bound) << endl;
      const double andrew_upper_bound = distinct_counts +
	upperbound_librarysize(andrew_counts_hist, max_terms);
      stats_out << "ANDREW_CF_UPPER=" << upper_bound << endl;
      stats_out << "ANDREW_CF_LOWER="
		<< andrew_cfa.lowerbound_librarysize(andrew_counts_hist, upper_bound) << endl;
      const double ztp_upper_bound = distinct_counts +
	upperbound_librarysize(ztp_counts_hist, max_terms);
      stats_out << "ZTP_CF_UPPER=" << upper_bound << endl;
      stats_out << "ZTP_CF_LOWER="
		<< ztp_cfa.lowerbound_librarysize(ztp_counts_hist, upper_bound) << endl;

		} */
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
