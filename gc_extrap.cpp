/*    gc_extrap: extrapolate genomic complexity 
 *
 *    Copyright (C) 2012 University of Southern California and
 *			 Andrew D. Smith and Timothy Daley
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
 *    along with this program.	If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <MappedRead.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>

#include "continued_fraction.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;


/*
 * This code is used to deal with read data in BAM format.
 */
/* dealing with bam version later
#ifdef HAVE_BAMTOOLS
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;
#endif
*/


/**************** FOR CLARITY BELOW WHEN COMPARING MAPPED READS **************/
static inline bool
same_strand(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() == b.get_strand();
}
static inline bool
strand_less(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() <= b.get_strand();
}
static inline bool
start_leq(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_start() <= b.get_start();
}
static inline bool
same_end(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_end() == b.get_end();
}
static inline bool
start_less(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_start() < b.get_start();
}
/******************************************************************************/


struct GenomicRegionOrderChecker {
  bool operator()(const GenomicRegion &prev, const GenomicRegion &gr) const {
    return end_two_check(prev, gr);
  }
  static bool 
  is_ready(const priority_queue<GenomicRegion, vector<GenomicRegion>, GenomicRegionOrderChecker> &pq,
	   const GenomicRegion &gr, const size_t max_width) {
    return !pq.top().same_chrom(gr) || pq.top().get_end() + max_width < gr.get_start();
  }
  static bool 
  end_two_check(const GenomicRegion &prev, const GenomicRegion &gr) {
    return (start_less(prev, gr) || 
	    (same_end(prev, gr) &&
	     (strand_less(prev, gr) ||
	      (same_strand(prev, gr) && start_leq(prev, gr)))));
  }
};



// add inputGR to the priority queue and if the priority queue is ready
static GenomicRegion
reorderGenomicRegions(const GenomicRegion &inputGR,
		      const size_t max_width,
		      priority_queue<GenomicRegion, vector<GenomicRegion>, 
				     GenomicRegionOrderChecker> &PQ){
  GenomicRegion outputGR;
  if(!PQ.empty() && GenomicRegionOrderChecker::is_ready(PQ, inputGR, max_width)){
    outputGR = PQ.top();
    PQ.pop();
  }

  PQ.push_back(inputGR);

  return outputGR;
}




// probabilistically split genomic regions into mutiple
// genomic regions of width equal to bin_size
static void
SplitGenomicRegion(const GenomicRegion &inputGR,
		       const size_t bin_size,
		       vector<GenomicRegion> &outputGRs){
  outputGRs.clear();
  GenomicRegion gr(inputGR);

  double frac =
    static_cast<double>(gr.get_start() % bin_size)/bin_size;
  const size_t width = gr.get_width();

  if(runif.runif(0.0, 1.0) > frac){
    gr.set_start(std::floor(static_cast<double>(gr.get_start())/
			    bin_size)*bin_size);
    gr.set_end(gr.get_start() + width);
  }
  else {
    gr.set_start(std::ceil(static_cast<double>(gr.get_start())/
			   bin_size)*bin_size);
    gr.set_end(gr.get_start() + width);
  }

  for(size_t i = 0; i < gr.get_width(); i += bin_size){

    const size_t curr_start = gr.get_start() + i;
    const size_t curr_end = std::min(gr.get_end(), curr_start + bin_size);
    frac = static_cast<double>(curr_end - curr_start)/bin_size;

    if(runif.runif(0.0, 1.0) <= frac){
      GenomicRegion binned_gr(gr.get_chrom(), curr_start, curr_end,
			      gr.get_name(), gr.get_score(), 
			      gr.get_strand());

      outputGRs.push_back(binned_gr);
    }
  }
}

/*
static size_t
load_values_MR_pe(const string input_file_name, vector<double> &values) {

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 MappedRead r, prev;
 if (!(in >> prev))
   throw "problem reading from: " + input_file_name;

 size_t n_reads = 1;
 values.push_back(1.0);
 while (in >> r) {
    if (r.same_chrom(prev) && r.get_start() < prev.get_start() 
	&& r.get_end() < prev.get_end())
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start() 
	|| r.get_end() != prev.get_end())
     values.push_back(1.0);
   else values.back()++;
   ++n_reads;
   prev.swap(r);
 }
 return n_reads;
}
*/

static size_t
load_values_MR_se(const string input_file_name, 
		  const size_t bin_size,
		  const size_t max_width,
		  vector<double> &values) {

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 MappedRead mr;
 if (!(in >> mr))
   throw "problem reading from: " + input_file_name;

 // initialize prioirty queue to reorder the split reads
 std::priority_queue<GenomicRegion, vector<GenomicRegion>, GenomicRegionOrderChecker> PQ;

 // prev and current Genomic Regions to compare
 GenomicRegion curr_gr, prev_gr;
 size_t n_reads = 1;
 values.push_back(1.0);
 do {
   vector<GenomicRegion> splitGRs;
   SplitGenomicRegion(mr.r,  bin_size, splitGRs);
   // add split Genomic Regions to the priority queue
   for(size_t i = 0; i < splitGRs.size(); i++)
     PQ.push(splitGRs[i]);

   // remove Genomic Regions from the priority queue
   while(!PQ.empty() && GenomicRegionOrderChecker::is_ready(PQ, splitGRs.back())){
     curr_gr = PQ.top();
     // only compare if the previous is not null (in the 1st iteration)
     if(prev_gr.get_chrom() != "(NULL)"){
       if(curr_gr.same_chrom(prev_gr) && curr_gr.get_start() > prev_gr.get_start()){
	 cerr << "current:\t" << curr_gr << endl;
	 cerr << "previous:\t" << prev_gr << endl;
	 throw SMITHLABException("split reads unsorted");
       }
       if(!curr_gr.same_chrom(prev_gr) || curr_gr.get_start() != prev.get_start())
	 values.push_back(1.0);
       else 
	 values.back()++;
     }
     prev_gr.swap(curr_gr);
   }

   n_reads++;
 } while (in >> mr);
 return n_reads;
}



void
resample_hist(const gsl_rng *rng, const vector<double> &vals_hist,
	      const double total_sampled_reads,
	      double expected_sample_size,
	      vector<double> &sample_hist) {
  
  const size_t hist_size = vals_hist.size();
  const double vals_mean = total_sampled_reads/expected_sample_size;
  
  sample_hist = vector<double>(hist_size, 0.0);
  vector<unsigned int> curr_sample(hist_size);
  double remaining = total_sampled_reads;
  
  while (remaining > 0) {
    
    // get a new sample
    expected_sample_size = max(1.0, (remaining/vals_mean)/2.0);
    gsl_ran_multinomial(rng, hist_size, 
			static_cast<unsigned int>(expected_sample_size),
			&vals_hist.front(), &curr_sample.front());
    
    // see how much we got
    double inc = 0.0;
    for (size_t i = 0; i < hist_size; ++i)
      inc += i*curr_sample[i];
    
    // only add to histogram if sampled reads < remaining reads
    if (inc <= remaining) {
      for (size_t i = 0; i < hist_size; i++)
	sample_hist[i] += static_cast<double>(curr_sample[i]);
      // update the amount we still need to get
      remaining -= inc;
    }
  }
}

static double
sample_count_distinct(const gsl_rng *rng,
		      const vector<size_t> &full_umis,
		      const size_t sample_size) {
  vector<size_t> sample_umis(sample_size);
  gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
		 (size_t *)&full_umis.front(), full_umis.size(), 
		 sizeof(size_t));
  double count = 1.0;
  for (size_t i = 1; i < sample_umis.size(); i++)
    if(sample_umis[i] != sample_umis[i-1])
      count++;

  return count;
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

void
estimates_bootstrap(const bool VERBOSE, const vector<double> &orig_values, 
		    const size_t bootstraps, const size_t orig_max_terms, 
		    const int diagonal, const double step_size, 
		    const double max_extrapolation, const double 
		    const double tolerance, const size_t max_iter,
		    vector<double> &Y_estimates,
		    vector< vector<double> > &yield_estimates) {
  // clear returning vectors
  yield_estimates.clear();
  
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
  const double max_val = max_extrapolation/vals_sum;
  
  for (size_t iter = 0; 
       (iter < 4*bootstraps && yield_estimates.size() < bootstraps); ++iter) {
    
    vector<double> yield_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist, vals_sum,
		  static_cast<double>(orig_values.size()), hist);

    const double initial_distinct = accumulate(hist.begin(), hist.end(), 0.0);
    
    //resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    //construct umi vector to sample from
    vector<size_t> umis;
    size_t umi = 1;
    for(size_t i = 1; i < hist.size(); i++){
      for(size_t j = 0; j < hist[i]; j++){
	for(size_t k = 0; k < i; k++)
	  umis.push_back(umi);
	umi++;
      }
    }
    assert(umis.size() == static_cast<size_t>(vals_sum));

    // compute complexity curve by random sampling w/out replacement
    size_t upper_limit = static_cast<size_t>(vals_sum);
    size_t step = static_cast<size_t>(step_size);
    size_t sample = step;
    while(sample < upper_limit){
      yield_vector.push_back(sample_count_distinct(rng, umis, sample));
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
    max_terms = max_terms - (max_terms % 2 == 0);
    
    //refit curve for lower bound
    const ContinuedFractionApproximation 
      lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
    
    const ContinuedFraction 
      lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));
    
    //extrapolate the curve start
    if (lower_cf.is_valid()){
      double sample_size = static_cast<double>(sample);
      while(sample_size < max_extrapolation){
	double t = (sample_size - vals_sum)/vals_sum;
	assert(t >= 0.0);
	yield_vector.push_back(initial_distinct + t*lower_cf(t));
	sample_size += step_size;
      }
    
    // SANITY CHECK    
      if (check_yield_estimates(yield_vector)) {
	yield_estimates.push_back(yield_vector);
	if (VERBOSE) cerr << '.';
	Y50_estimates.push_back(lower_cf.Y50(hist, vals_sum, max_val, tolerance, max_iter));
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
  if (yield_estimates.size() < bootstraps)
    throw SMITHLABException("too many iterations, poor sample");
}


static inline double
alpha_log_confint_multiplier(const double estimate,
			     const double variance, const double alpha) {
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  return exp(inv_norm_alpha*
	     sqrt(log(1.0 + variance/pow(estimate, 2))));
}


static void
vector_median_and_ci(const vector<vector<double> > &estimates,
		     const double alpha, 
		     vector<double> &median_estimates,
		     vector<double> &lower_ci_lognormal, 
		     vector<double> &upper_ci_lognormal) {
  
  assert(!estimates.empty());
  
  const size_t n_est = estimates.size();
  vector<double> estimates_row(estimates.size(), 0.0);
  for (size_t i = 0; i < estimates[0].size(); i++) {
    
    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = estimates[k][i];
    
    // sort to get confidence interval
    sort(estimates_row.begin(), estimates_row.end());
    const double curr_median = 
      gsl_stats_median_from_sorted_data(&estimates_row[0], 1, n_est);
    
    median_estimates.push_back(curr_median);
    const double variance = gsl_stats_variance(&estimates_row[0], 1, n_est);
    const double confint_mltr = 
      alpha_log_confint_multiplier(curr_median, variance, alpha);
    lower_ci_lognormal.push_back(curr_median/confint_mltr);
    upper_ci_lognormal.push_back(curr_median*confint_mltr);
  }
}

static void
median_and_ci(const vector<double> &estimates,
	      const double alpha,
	      double &median_estimate,
	      double &lower_ci_estimate,
	      double &upper_ci_estimate){
  assert(!estimates.empty());
  const size_t n_est = estimates.size();
  vector<double> sorted_estimates(estimates);
  sort(sorted_estimates.begin(), sorted_estimates.end());
  median_estimate = 
    gsl_stats_median_from_sorted_data(&sorted_estimates[0], 1, n_est);
  const double variance = gsl_stats_variance(&sorted_estimates[0], 1, n_est);
  const double confint_mltr = 
    alpha_log_confint_multiplier(median_estimate, variance, alpha);

  lower_ci_estimate = median_estimate/confint_mltr;
  upper_ci_estimate = median_estimate*confint_mltr;

}

static void
write_predicted_curve(const string outfile, const double values_sum,
		      const double c_level, const double step_size,
		      const size_t bin_size,
		      const vector<double> &median_yield_estimates,
		      const vector<double> &yield_lower_ci_lognormal,
		      const vector<double> &yield_upper_ci_lognormal) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out << "TOTAL_READS\tEXPECTED_COVERED_BASES\t"
      << "LOWER_" << 100*c_level << "%CI\t"
      << "UPPER_" << 100*c_level << "%CI" << endl;
  
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);
  
  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < median_yield_estimates.size(); ++i)
    out << (i + 1)*step_size*bin_size << '\t' 
	<< median_yield_estimates[i]*bin_size << '\t'
	<< yield_lower_ci_lognormal[i]*bin_size << '\t' 
	<< yield_upper_ci_lognormal[i]*bin_size << endl;
}



int
main(const int argc, const char **argv) {

  try {
    
    const size_t MIN_REQUIRED_COUNTS = 8;

    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 1000;
    double max_extrapolation = 1.0e10;
    double step_size = 1e7;
    size_t bootstraps = 100;
    int diagonal = -1;
    double c_level = 0.95;
    double tolerance = 1e-20;
    size_t max_width = 1000;
    size_t max_iter = 100;
    double tolerance = 1.0e-20;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool PAIRED_END = false;
    
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
		      false , outfile);
    opt_parse.add_opt("max_width", 'w', "max fragment length, "
		      "set equal to read length for single end reads",
		      false, max_width);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation "
		      "(default: " + toa(max_extrapolation) + ")",
		      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size in bases between extrapolations "
		      "(default: " + toa(step_size) + ")", 
		      false, step_size);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps "
		      "(default: " + toa(bootstraps) + "), ",
		      false, bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals "
		      "(default: " + toa(c_level) + ")", false, c_level);
    opt_parse.add_opt("terms",'x',"maximum number of terms", 
		      false, orig_max_terms);
    //    opt_parse.add_opt("tol",'t', "numerical tolerance", false, tolerance);
    //    opt_parse.add_opt("max_iter",'i', "maximum number of iteration",
    //		      false, max_iter);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false, VERBOSE);
    //    opt_parse.add_opt("pe", 'P', "input is paired end read file",
    //		      false, PAIRED_END);
    
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
    size_t n_reads = 0;

    /*    if(PAIRED_END)
      load_values_BED_pe(input_file_name, values);  
      else */
      n_reads = load_values_BED_se(input_file_name, values);

    
    // JUST A SANITY CHECK
    const double total_reads = accumulate(values.begin(), values.end(), 0.0);

    // for large initial experiments need to adjust step size
    // otherwise small relative steps do not account for variance
    // in extrapolation
    if(step_size < (total_reads/20)){
       step_size = std::max(1e6, 1e6*round(total_reads/(20*step_size)));
       if(VERBOSE)
	 cerr << "ADJUSTED_STEP_SIZE = " << step_size << endl;
    }
       
    const size_t max_observed_count = 
      static_cast<size_t>(*std::max_element(values.begin(), values.end()));


    
    // BUILD THE HISTOGRAM
    vector<double> counts_hist(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_hist[static_cast<size_t>(values[i])];
    
    const size_t distinct_counts = 
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
					bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
      cerr << "TOTAL READS         = " << n_reads << endl
	   << "DISTINCT BINS       = " << total_reads << endl
	   << "TOTAL BASES         = " << total_reads*bin_size << endl
	   << "TOTAL COVERED BASES = " << values.size()*bin_size << endl
	   << "MAX COUNT           = " << max_observed_count << endl
	   << "COUNTS OF 1         = " << counts_hist[1] << endl;
    
    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    }

    // catch if all reads are distinct
    if (max_observed_count < MIN_REQUIRED_COUNTS)
      throw SMITHLABException("sample not sufficiently deep or duplicates removed");
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // BOOTSTRAPS

    if(bootstraps < 10)
      throw SMITHLABException("too few bootstraps, must be at least 10");

    if (VERBOSE) 
      cerr << "[BOOTSTRAP ESTIMATES]" << endl;
      
    vector<vector <double> > yield_estimates;
    //    vector< vector<double> > sat_estimates;
    vector<double> Y50_estimates;
    vector<double> lower_libsize, upper_libsize;
    estimates_bootstrap(VERBOSE, values,  bootstraps, orig_max_terms,
			diagonal, step_size, max_extrapolation, tolerance,
			max_iter, yield_estimates);
      
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    if (VERBOSE)
      cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
 
    // yield vector median and ci    
    vector<double> median_yield_estimates;
    vector<double> yield_upper_ci_lognormal, yield_lower_ci_lognormal,
      yield_upper_ci_quantile, yield_lower_ci_quantile;
    vector_median_and_ci(yield_estimates, 1.0 - c_level, 
			 median_yield_estimates, 
			 yield_lower_ci_lognormal, yield_upper_ci_lognormal);



    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    if (VERBOSE) 
      cerr << "[WRITING OUTPUT]" << endl;
    
    write_predicted_curve(outfile, total_reads, c_level, step_size,
			  median_yield_estimates,
			  yield_lower_ci_lognormal, yield_upper_ci_lognormal);


      
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
