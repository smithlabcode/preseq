/*    saturation_extrap:
 *
 *    Copyright (C) 2013 University of Southern California and
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

static GenomicRegion
BamToGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
		   const BamAlignment &ba){

  const unordered_map<size_t, string>::const_iterator
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = ba.Position + ba.InsertSize;

  return GenomicRegion(chrom, start, end);

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
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
     //only count unpaired reads or the left mate of paired reads
      if(!(bam.IsPaired()) || 
	 (bam.IsFirstMate())){

	SimpleGenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam));
	if (r.same_chrom(prev) && r.get_start() < prev.get_start())
	  throw SMITHLABException("locations unsorted in: " + input_file_name);
    
	if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
	  values.push_back(1.0);
	else values.back()++;
	++n_reads;
	prev.swap(r);
      }
    }
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

  GenomicRegion prev;
  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
      // ignore reads that do not map concoordantly
      if(bam.IsPaired() && bam.IsProperPair() && bam.IsFirstMate()){
	GenomicRegion r(BamToGenomicRegion(chrom_lookup, bam));
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
    }
  }
  reader.Close();
  return n_reads;
}
#endif


static size_t
load_values_BED_se(const string input_file_name, vector<double> &values) {

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 SimpleGenomicRegion r, prev;
 if (!(in >> prev))
   throw "problem reading from: " + input_file_name;

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
    gsl_ran_multinomial(rng, hist_size, expected_sample_size,
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
bootstrap_saturation_deriv(const bool VERBOSE, const vector<double> &orig_values,
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

    vector<double> hist;

    resample_hist(rng, orig_hist, vals_sum,
		  static_cast<double>(orig_values.size()),
		  hist);
    
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
      lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));

    vector<double> saturation_estimates;
    if(lower_cf.is_valid()){
      lower_cf.extrapolate_yield_deriv(hist, vals_sum,
				      max_val, val_step,
				      saturation_estimates);

    }
    else if(VERBOSE)
      cerr << "not_valid" << endl;

    if(check_saturation_estimates(saturation_estimates)){
      full_estimates.push_back(saturation_estimates);
      if(VERBOSE)
	cerr << ".";
    }
    else if(VERBOSE)
      cerr << "fail_check" << endl;

  }
  cerr << "deriv estimates size =" << full_estimates.size();
}

/*
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
    resample_hist(rng, orig_hist, vals_sum,
		  static_cast<double>(orig_values.size()),
		  hist);
    
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
    max_terms = max_terms - (max_terms % 2 == 1);
    
    //refit curve for lower bound
    const ContinuedFractionApproximation 
      lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
    
    const ContinuedFraction 
      lower_cf(lower_cfa.optimal_cont_frac_satur(hist));

    vector<double> saturation_estimates;
    if(lower_cf.is_valid()){
      lower_cf.extrapolate_saturation(hist, vals_sum,
				      max_val, val_step,
				      saturation_estimates);
      saturation_estimates.insert(saturation_estimates.begin(),
				  hist[1]/vals_sum);
    }
    else if(VERBOSE)
      cerr << "not_valid" << endl;
    if(check_saturation_estimates(saturation_estimates)){
      full_estimates.push_back(saturation_estimates);
      if(VERBOSE) 
	cerr << ".";
    }
    else if(VERBOSE)
      cerr << "fail_check" << endl;

  }

  cerr << "cf estimates size = " << full_estimates.size() << endl;
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
    lower_cf(lower_cfa.optimal_cont_frac_yield(hist));
  if(lower_cf.is_valid())
    lower_cf.extrapolate_saturation(hist, vals_sum,
				    max_val, val_step,
				    estimates);
  if(!check_saturation_estimates(estimates))
    cerr << "ESTIMATES UNSTABLE, MORE DATA REQUIRED" << endl; 

}
*/

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
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    
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
    
    // READ IN THE DATA
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

    const double vals_sum = accumulate(values.begin(), values.end(), 0.0);
    
    const double max_val = max_extrapolation/vals_sum;
    const double val_step = step_size/vals_sum;
    
    const size_t max_observed_count = 
      static_cast<size_t>(*std::max_element(values.begin(), values.end()));

    // catch if all reads are distinct
    if(max_observed_count < 8)
      throw SMITHLABException("sample not sufficiently distinct, unable to estimate");
    
    // BUILD THE HISTOGRAM
    vector<double> counts_hist(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_hist[static_cast<size_t>(values[i])];
    
    if (VERBOSE)
      cerr << "TOTAL READS     = " << vals_sum << endl
	   << "DISTINCT READS  = " << values.size() << endl
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

      vector< vector<double> > full_deriv_estimates;
      bootstrap_saturation_deriv(VERBOSE, values,  
			   bootstraps, orig_max_terms, diagonal, 
			   step_size, max_extrapolation, 
			   max_val, val_step, 
			   full_deriv_estimates);



      if(VERBOSE)
	cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;

      vector<double> mean_deriv_estimates, upper_alphaCI_deriv, lower_alphaCI_deriv;
      compute_mean_and_alphaCI(full_deriv_estimates, 1.0 - c_level,
			       mean_deriv_estimates, lower_alphaCI_deriv,
			       upper_alphaCI_deriv);


      if(VERBOSE)
	cerr << "[WRITING OUTPUT]" << endl;

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      out << "TOTAL_READS" << '\t' 
	  << "SATURATION" << '\t' 
	  << "LOWER_" << 100*c_level << "%CI" << '\t' 
	  << "UPPER_" << 100*c_level << "%CI" << endl;
    
      double val = 0.0;
      for (size_t i = 0; i < mean_deriv_estimates.size(); ++i, val += val_step)
	out << (val + 1.0)*vals_sum << '\t' 
	    << mean_deriv_estimates[i] << '\t'
	    << lower_alphaCI_deriv[i] << '\t' 
	    << upper_alphaCI_deriv[i] << endl;
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
