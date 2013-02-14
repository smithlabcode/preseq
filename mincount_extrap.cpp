/*    lc_extrap:
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

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
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


static size_t
load_values_BAM_se(const string &input_file_name, vector<double> &values) {
  
  BamReader reader;
  reader.Open(input_file_name);

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;

  size_t n_reads = 1;
  values.push_back(1.0);

  SimpleGenomicRegion prev;
  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    SimpleGenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam));
    if (r.same_chrom(prev) && r.get_start() < prev.get_start())
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
      values.push_back(1.0);
    else values.back()++;
    ++n_reads;
    prev.swap(r);
  }
  reader.Close();
  return n_reads;
}

static size_t
load_values_BAM_pe(const string &input_file_name, vector<double> &values) {

  BamReader reader;
  reader.Open(input_file_name);

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;

  size_t n_reads = 1;
  values.push_back(1.0);

  SimpleGenomicRegion prev;
  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    SimpleGenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam));
    if (r.same_chrom(prev) && r.get_start() < prev.get_start())
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
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
load_values_BED_se(const string input_file_name, vector<double> &values) {
  
  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " + input_file_name);
  
  SimpleGenomicRegion r, prev;
  if (!(in >> prev))
    throw SMITHLABException("problem reading from: " + input_file_name);
  
  size_t n_reads = 1;
  values.push_back(1.0);
  while (in >> r) {
    if (r.same_chrom(prev) && r.get_start() < prev.get_start())
      throw SMITHLABException("locations unsorted in: " + input_file_name);

    if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
      values.push_back(1.0);
    else values.back()++;
    ++n_reads;
    prev.swap(r);
  }
  return n_reads;
}

static size_t
load_values_BED_pe(const string input_file_name, vector<double> &values) {

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 GenomicRegion r, prev;
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


static size_t
load_values(const string input_file_name, vector<double> &values) {

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " + input_file_name);

  vector<double> full_values;
  size_t n_reads = 0;
  static const size_t buffer_size = 10000; // Magic!
  while(!in.eof()){
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    full_values.push_back(atof(buffer));
    if(full_values.back() <= 0.0){
      cerr << "INVALID INPUT\t" << buffer << endl;
      throw SMITHLABException("ERROR IN INPUT");
    }
    ++n_reads;
    in.peek();
  }
  in.close();
  if(full_values.back() == 0)
    full_values.pop_back();

  values.swap(full_values);
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
sample_count_reads_w_mincount(const gsl_rng *rng,
			      const vector<size_t> &full_umis,
			      const size_t sample_size,
			      const size_t mincount) {
  vector<size_t> sample_umis(sample_size);
  gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
		 (size_t *)&full_umis.front(), full_umis.size(), 
		 sizeof(size_t));

  double number_observed = 0.0;
  size_t current_count = 1;

  for (size_t i = 1; i < sample_umis.size(); i++){
    if(sample_umis[i] == sample_umis[i-1])
      current_count++;
    else{
      if(current_count >= mincount)
	number_observed++;
      current_count = 1;
    }
  }
  if(current_count >= mincount)
    number_observed++;

  return number_observed;
}


static bool
check_mincount_estimates_stability(const vector<double> &estimates,
				   const double max_change_per_time_step) {
  // make sure that the estimate is increasing in the time_step and
  // is below the initial distinct per step_size
  for (size_t i = 1; i < estimates.size(); ++i){
    if(!finite(estimates[i])){
      return false;
    }
    if ((estimates[i] < estimates[i - 1]) ||
	(estimates[i] - estimates[i - 1] > max_change_per_time_step)){
      return false;
    }
  }    

  return true;
}


void
estimates_bootstrap(const bool VERBOSE, const vector<double> &orig_values, 
		    const size_t bootstraps, const size_t orig_max_terms, 
		    const int order, const double step_size, 
		    const double max_extrapolation, const size_t mincount,
		    vector< vector<double> > &full_estimates) {
  // clear returning vectors
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
  
  const double vals_sum = accumulate(orig_values.begin(), orig_values.end(), 0.0);
  
  for (size_t iter = 0; 
       (iter < 2*bootstraps && full_estimates.size() < bootstraps); ++iter) {
    
    vector<double> mincount_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist, vals_sum,
		  static_cast<double>(orig_values.size()), hist);

    /*   cerr << "sampled hist:" << endl;
    for(size_t i = 0; i <= std::min(hist.size(), orig_max_terms); i++)
      if(hist[i] > 0)
	cerr << i << "\t" << hist[i] << endl;
    */

    const double initial_observed = accumulate(hist.begin() + mincount, hist.end(), 0.0);
    
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
      mincount_vector.push_back(sample_count_reads_w_mincount(rng, umis, sample, mincount));
      sample += step;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < hist.size() &&
	   hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;   

    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - mincount);
    max_terms = max_terms - (max_terms % 2 == 0);

    const ContinuedFractionApproximation 
      mincount_cfa(order, max_terms, step_size, max_extrapolation);
    
    const ContinuedFraction 
      mincount_cf(mincount_cfa.optimal_cont_frac_mincount(hist, mincount, order));


    //extrapolate the curve start
    if (mincount_cf.is_valid()){
      //  cerr << "valid" << endl;
      double sample_size = static_cast<double>(sample);
      while(sample_size <= max_extrapolation){
	double t = (sample_size - vals_sum)/vals_sum;
	assert(t >= 0.0);
	mincount_vector.push_back(initial_observed + t*mincount_cf(t));
	sample_size += step_size;
      }
      if(check_mincount_estimates_stability(mincount_vector, step_size)){
	full_estimates.push_back(mincount_vector);
	if (VERBOSE) cerr << '.';
      }
      else if(VERBOSE) cerr << '_';
    }
  else if(VERBOSE) cerr << '_';
  }
  if (VERBOSE)
    cerr << endl;
  if (full_estimates.size() < bootstraps)
    throw SMITHLABException("too many iterations, poor sample");
}


static void
return_median_and_ci(const vector<vector<double> > &estimates,
		     const double alpha, const double initial_distinct,
		     const double vals_sum, const double step_size, 
		     vector<double> &median_estimates,
		     vector<double> &lower_ci, vector<double> &upper_ci) {
  //  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
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
    // const double std_dev = sqrt(gsl_stats_variance(&estimates_row[0], 1, n_est));
    upper_ci.push_back(gsl_stats_quantile_from_sorted_data(&estimates_row[0], 1, n_est, 0.95));
    lower_ci.push_back(gsl_stats_quantile_from_sorted_data(&estimates_row[0], 1, n_est, 0.05));
  }
}


static void
write_predicted_curve(const string outfile, const double c_level,
		      const double step_size,
		      const vector<double> &median_yield_estimates,
		      const vector<double> &yield_lower_ci,
		      const vector<double> &yield_upper_ci) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOWER_" << 100*c_level << "%CI\t"
      << "UPPER_" << 100*c_level << "%CI" << endl;
  
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);
  
  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < median_yield_estimates.size(); ++i)
    out << (i + 1)*step_size << '\t' 
	<< median_yield_estimates[i] << '\t'
	<< yield_lower_ci[i] << '\t' << yield_upper_ci[i] << endl;
}



int
main(const int argc, const char **argv) {

  try {
    
    const size_t MIN_REQUIRED_COUNTS = 8;

    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 200;
    double max_extrapolation = 1.0e10;
    double step_size = 1e6;
    size_t bootstraps = 100;
    int order = 0;
    double c_level = 0.95;
    size_t mincount = 2;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    
#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
		      false , outfile);
    opt_parse.add_opt("extrap",'e',"maximum extrapolation "
		      "(default: " + toa(max_extrapolation) + ")",
		      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size in extrapolations "
		      "(default: " + toa(step_size) + ")", 
		      false, step_size);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps "
		      "(default: " + toa(bootstraps) + "), ",
		      false, bootstraps);
    opt_parse.add_opt("cval", 'c', "level for confidence intervals "
		      "(default: " + toa(c_level) + ")", false, c_level);
    opt_parse.add_opt("mincount",'m',"minimum number of observations for a read to be counted "
		      "(default: " + toa(mincount) + "), ",
		      false, mincount); 
    //	opt_parse.add_opt("terms",'t',"maximum number of terms", false, 
    //	     orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false, VERBOSE);
#ifdef HAVE_BAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format", 
		      false, BAM_FORMAT_INPUT);
#endif
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
		      false, PAIRED_END);
    opt_parse.add_opt("vals", 'V', 
		      "input is a text file containing only the observed counts",
		      false, VALS_INPUT);
    
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
    if(VALS_INPUT)
      load_values(input_file_name, values);
#ifdef HAVE_BAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END)
      load_values_BAM_pe(input_file_name, values);
    else if(BAM_FORMAT_INPUT)
      load_values_BAM_se(input_file_name, values);
#endif
    else if(PAIRED_END)
      load_values_BED_pe(input_file_name, values);  
    else
      load_values_BED_se(input_file_name, values);
    
    // JUST A SANITY CHECK
    const double values_sum = accumulate(values.begin(), values.end(), 0.0);
        
    const size_t max_observed_count = 
      static_cast<size_t>(*std::max_element(values.begin(), values.end()));

    // catch if all reads are distinct
    if (max_observed_count < MIN_REQUIRED_COUNTS)
      throw SMITHLABException("sample appears too uniform");
    
    // BUILD THE HISTOGRAM
    vector<double> counts_hist(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_hist[static_cast<size_t>(values[i])];

    const double initial_observed =
      accumulate(counts_hist.begin() + mincount, counts_hist.end(), 0.0);

    if (VERBOSE)
      cerr << "TOTAL READS      = " << values_sum << endl
	   << "DISTINCT READS   = " << values.size() << endl
	   << "INITIAL OBSERVED = " << initial_observed << endl
	   << "MAX COUNT        = " << max_observed_count << endl
	   << "COUNTS OF 1      = " << counts_hist[1] << endl
	   << "MAX TERMS        = " << orig_max_terms << endl;
    
    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      /*
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
      */
    }
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // BOOTSTRAPS

    //    if(bootstraps < 10)
    //  throw SMITHLABException("too few bootstraps, must be at least 10");

    if (VERBOSE) 
      cerr << "[BOOTSTRAP ESTIMATES]" << endl;
      
    vector<vector <double> > yield_estimates;
    vector< vector<double> > sat_estimates;
    vector<double> lower_libsize, upper_libsize;
    estimates_bootstrap(VERBOSE, values,  bootstraps, orig_max_terms,
			order, step_size, max_extrapolation, 
			mincount, yield_estimates);
      
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // BOOTSTRAPS
    if (VERBOSE)
      cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
      
    vector<double> median_yield_estimates;
    vector<double> yield_upper_ci, yield_lower_ci;
    return_median_and_ci(yield_estimates, 1.0 - c_level, 
			 initial_observed, values_sum, step_size,
			 median_yield_estimates, 
			 yield_lower_ci, yield_upper_ci);

    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    if (VERBOSE) 
      cerr << "[WRITING OUTPUT]" << endl;
    
    write_predicted_curve(outfile, c_level, step_size,
			  median_yield_estimates,
			  yield_lower_ci, yield_upper_ci);
      
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
