/*    saturation_extrap:
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

#include "continued_fraction.hpp"

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

static SimpleGenomicRegion
BamToSimpleGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
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
ReadBAMFormatInput(const string &infile, 
		   vector<SimpleGenomicRegion> &read_locations) {
  
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
    read_locations.push_back(BamToSimpleGenomicRegion(chrom_lookup, bam));
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


/*
void
resample_values(const vector<double> &in_values,
		const gsl_rng *rng,
		vector<double> &out_values) {
  out_values.clear();
  const double orig_number_captures = accumulate(in_values.begin(), 
						 in_values.end(), 0.0);
  double sampled_number_captures = 0.0;
  while(sampled_number_captures < orig_number_captures) {
    double u = gsl_rng_uniform(rng);
    while(u == 1.0)
      u = gsl_rng_uniform(rng);
    const size_t sample_indx  = 
      static_cast<size_t>(floor(in_values.size()*u));
    sampled_number_captures += in_values[sample_indx];
    if (sampled_number_captures > orig_number_captures)
      out_values.push_back(static_cast<double>(sampled_number_captures
					       - orig_number_captures));
    // ensure that number of captures is constant in resampling
    else
      out_values.push_back(in_values[sample_indx]);
  } 
}
*/

/*void
resample_values(const vector<double> &full_values,
		const gsl_rng *rng,
		const size_t estimated_sample_size,
		const double total_sampled_reads,
		vector<double> &sample_values){
  sample_values.clear();
  vector<double> temp_sample(estimated_sample_size, 0.0);
  gsl_ran_sample(rng, (double *)&temp_sample.front(), estimated_sample_size,
		 (double *)&full_values.front(), full_values.size(), sizeof(double));
  double temp_sample_sum = accumulate(temp_sample.begin(), temp_sample.end(), 0.0);

  // too few reads, add some
  if(temp_sample_sum < estimated_sample_size){
    while(temp_sample_sum < estimated_sample_size){
      temp_sample.push_back(full_values[rand() % full_values.size()]);
      temp_sample_sum += temp_sample.back();
    }
  }
  // too many reads, delete some
  if(temp_sample_sum > estimated_sample_size){
    while(temp_sample_sum > estimated_sample_size){
      double deleted_val = temp_sample.back();
      temp_sample.pop_back();
      temp_sample_sum -= deleted_val;
    }
  }
  
  // just right
  temp_sample.push_back(total_sampled_reads - temp_sample_sum);
  // assert(accumulate(temp_sample.begin(), temp_sample.end(), 0.0) == total_sampled_reads);
  sample_values.swap(temp_sample);
}
*/

//sample histogram
void
resample_hist(const vector<double> &values,
	      const vector<double> &vals_hist,
	      const gsl_rng *rng,
	      const double total_sampled_reads,
	      const size_t expected_sample_size,
	      vector<double> &sample_hist){
  double proportions[vals_hist.size()];
  for(size_t i = 0; i < vals_hist.size(); i++)
    proportions[i] = vals_hist[i];
  //erase the 0th entry
  unsigned int temp_hist[vals_hist.size()];
  gsl_ran_multinomial(rng, vals_hist.size(), expected_sample_size,
		      proportions, temp_hist);

  sample_hist.clear();
  for(size_t i =0; i < vals_hist.size(); i++)
    sample_hist.push_back(static_cast<double>(temp_hist[i]));

  double sample_hist_sum = 0.0;
  for(size_t i = 0; i < sample_hist.size(); i++)
    sample_hist_sum += i*sample_hist[i];

  // too few reads, add some
  if(sample_hist_sum < total_sampled_reads){
    while(sample_hist_sum < total_sampled_reads){
      double temp_val = values[rand() % values.size()];
      sample_hist[static_cast<size_t>(temp_val)]++;
      sample_hist_sum += temp_val;
    }
  }

  // too many_reads, delete some
  if(sample_hist_sum > total_sampled_reads){
    while(sample_hist_sum > total_sampled_reads){
      size_t temp_val = static_cast<size_t>(values[rand() % values.size()]);
      if(sample_hist[temp_val] > 0){
	sample_hist[temp_val]--;
	sample_hist_sum -= static_cast<double>(temp_val);
      }
    }
  }

  // just right
  assert(sample_hist_sum <= total_sampled_reads);
  sample_hist[static_cast<size_t>(total_sampled_reads - sample_hist_sum)]++;
  sample_hist_sum += total_sampled_reads - sample_hist_sum;

  assert(sample_hist_sum == total_sampled_reads);
}


static inline bool
check_saturation_estimates(const vector<double> estimates){
  if(estimates.empty())
    return false;

  // make sure estimates are decreasing and
  // between 0 & 1
  if(estimates[0] >= 1.0 || estimates[0] < 0.0)
    return false;

  for(size_t i = 1; i < estimates.size(); i++)
    if(estimates[i] > estimates[i-1] ||
       estimates[i] >= 1.0 ||
       estimates[i] < 0.0) 
      return false;
  
  return true;
}

void
bootstrap_saturation(const bool VERBOSE, const vector<double> &orig_values,
		     const size_t bootstraps, const size_t orig_max_terms, 
		     const int diagonal, const double step_size, 
		     const double max_extrapolation, 
		     const double max_val, const double val_step, 
		     vector< vector<double> > &full_estimates){
  full_estimates.clear();

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


  const double vals_sum =
     accumulate(orig_values.begin(), orig_values.end(), 0.0);

  for (size_t iter = 0; 
       iter < 2*bootstraps && full_estimates.size() < bootstraps; ++iter) {
    
    //    vector<double> boot_values;
    //  resample_values(orig_values, rng, orig_values.size(), vals_sum, boot_values);
    
    //  const size_t max_observed_count = 
    //     static_cast<size_t>(*std::max_element(boot_values.begin(), 
    //					    boot_values.end()));
    
    // BUILD THE HISTOGRAM
    //  vector<double> hist(max_observed_count + 1, 0.0);
    //  for (size_t i = 0; i < boot_values.size(); ++i)
    //   ++hist[static_cast<size_t>(boot_values[i])];

    vector<double> hist;
    resample_hist(orig_values, orig_hist, rng, vals_sum, 
		  orig_values.size(), hist);
    
    //resize boot_hist to remove excess zeros
    while(hist.back() == 0)
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
      lower_cf(lower_cfa.optimal_continued_fraction(hist));

    vector<double> saturation_estimates;
    if(lower_cf.is_valid())
      lower_cf.extrapolate_saturation(hist, vals_sum,
				      max_val, val_step,
				      saturation_estimates);
    if(check_saturation_estimates(saturation_estimates))
      full_estimates.push_back(saturation_estimates);

  }
}

void 
single_estimates(const bool VERBOSE, const vector<double> &hist,
		 const double vals_sum,
		 size_t max_terms, const int diagonal,
		 const double step_size, const double max_extrapolation,
		 const double max_val, const double val_step, 
		 vector<double> &estimates){
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
    lower_cf(lower_cfa.optimal_continued_fraction(hist));
  if(lower_cf.is_valid())
    lower_cf.extrapolate_saturation(hist, vals_sum,
				    max_val, val_step,
				    estimates);
  if(!check_saturation_estimates(estimates))
    cerr << "ESTIMATES UNSTABLE, MORE DATA REQUIRED" << endl; 

}

static double
compute_var(const vector<double> &estimates,
	    const double mean) {
  double variance = 0.0;
  for(size_t i = 0; i < estimates.size(); i++)
    variance += (estimates[i] - mean)*(estimates[i] - mean)/estimates.size();
  return variance; 
}


static void
compute_mean_and_alphaCI(const vector< vector<double> > &full_estimates,
			 const double alpha, 
			 vector<double> &mean_estimates,
			 vector<double> &lower_alphaCI,
			 vector<double> &upper_alphaCI){
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);

  for(size_t i = 0; i < full_estimates[0].size(); i++) {
    // estimates is in wrong order, work locally on const val
    vector<double> log_estimates_row(full_estimates.size(), 0.0);

    for(size_t k = 0; k < log_estimates_row.size(); ++k)
      log_estimates_row[k] = log(full_estimates[k][i]);

    mean_estimates.push_back(exp(accumulate(log_estimates_row.begin(),
					    log_estimates_row.end(), 0.0)/
				 log_estimates_row.size()));
    const double variance = compute_var(log_estimates_row,
					log(mean_estimates.back()));

    // log confidence intervals
    upper_alphaCI.push_back(exp(log(mean_estimates.back())
				    + inv_norm_alpha*sqrt(variance)));
    lower_alphaCI.push_back(exp(log(mean_estimates.back())
				    - inv_norm_alpha*sqrt(variance)));

  }

}


int
main(const int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 200;
    double max_extrapolation = 1.0e10;
    double step_size = 1e6;
    size_t bootstraps = 100;
    int diagonal = -1;
    double c_level = 0.95;

    
    /* FLAGS */
    bool VERBOSE = false;
    
#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "",
			   "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
		      false , outfile);
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
    //    opt_parse.add_opt("terms",'t',"maximum number of terms", false, 
    //		      orig_max_terms);
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
    
    // JUST A SANITY CHECK
    const size_t n_reads = read_locations.size();
    const size_t values_sum = 
      static_cast<size_t>(accumulate(values.begin(), values.end(), 0.0));
    assert(values_sum == n_reads);
    
    const double max_val = max_extrapolation/static_cast<double>(values_sum);
    const double val_step = step_size/static_cast<double>(values_sum);
    
    const size_t max_observed_count = 
      static_cast<size_t>(*std::max_element(values.begin(), values.end()));

    // catch if all reads are distinct
    if(max_observed_count < 8)
      throw SMITHLABException("sample not sufficiently distinct, unable to estimate");
    
    // BUILD THE HISTOGRAM
    vector<double> counts_hist(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_hist[static_cast<size_t>(values[i])];
    
    const size_t distinct_counts = 
      static_cast<size_t>(count_if(counts_hist.begin(), counts_hist.end(),
				   bind2nd(std::greater<double>(), 0.0)));
    
    if (VERBOSE)
      cerr << "TOTAL READS     = " << read_locations.size() << endl
	   << "DISTINCT COUNTS = " << distinct_counts << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_hist[1] << endl
	   << "MAX TERMS       = " << orig_max_terms << endl;
    
    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for(size_t i = 0; i < counts_hist.size(); i++)
	if(counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    }

    bool BOOTSTRAP = (bootstraps > 10);
    if(BOOTSTRAP){
      if(VERBOSE)
	cerr << "[BOOTSTRAP ESTIMATES]" << endl;

      vector< vector<double> > full_estimates;
      bootstrap_saturation(VERBOSE, values,  
			   bootstraps, orig_max_terms, diagonal, 
			   step_size, max_extrapolation, 
			   max_val, val_step, 
			   full_estimates);

      if(VERBOSE)
	cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;

      vector<double> mean_estimates, upper_alphaCI, lower_alphaCI;
      compute_mean_and_alphaCI(full_estimates, 1.0 - c_level,
			       mean_estimates, lower_alphaCI,
			       upper_alphaCI);

      if(VERBOSE)
	cerr << "[WRITING OUTPUT]" << endl;

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "#TOTAL_READS" << '\t' 
	  << "EXPECTED_SATURATION" << '\t' 
	  << "LOWER_" << 100*c_level << "%CI" << '\t' 
	  << "UPPER_" << 100*c_level << "%CI" << endl;
    
      double val = 0.0;
      for (size_t i = 0; i < mean_estimates.size(); ++i, val += val_step)
	out << (val + 1.0)*values_sum << '\t' 
	    << mean_estimates[i] << '\t'
	    << lower_alphaCI[i] << '\t'
	    << upper_alphaCI[i] << endl;
    }
    else{
      if(VERBOSE)
	cerr << "[ESTIMATING]" << endl;
      vector<double> estimates;
      single_estimates(VERBOSE, counts_hist, values_sum, orig_max_terms, 
		       diagonal, step_size,  max_extrapolation, 
		       max_val, val_step, estimates);

      if (VERBOSE) 
	cerr << "[WRITING OUTPUT]" << endl;
    
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
      out << "#TOTAL_READS" << '\t' 
	  << "EXPECTED_SATURATION" <<  endl;
    
      double val = 0.0;
      for (size_t i = 0; i < estimates.size(); ++i, val += val_step)
	out << (val + 1.0)*values_sum << '\t' 
	    << estimates[i] << endl;

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
