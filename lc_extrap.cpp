/*    lc_extrap: extrapolate complexity curve
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
    double val = atof(buffer);
    if(val > 0.0)
      full_values.push_back(val);
    if(full_values.back() < 0.0){
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

static void
load_histogram(const string &filename, vector<double> &hist) {

  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("could not open histogram: " + filename);

  size_t line_count = 0ul, prev_read_count = 0ul;
  string buffer;
  while (getline(in, buffer)) {
    ++line_count;
    size_t read_count = 0ul;
    double frequency = 0.0;
    std::istringstream is(buffer);
    if (!(is >> read_count >> frequency))
      throw SMITHLABException("bad histogram line format:\n" +
                              buffer + "\n(line " + toa(line_count) + ")");
    if (read_count < prev_read_count)
      throw SMITHLABException("bad line order in file " +
                              filename + "\n(line " +
                              toa(line_count) + ")");
    hist.resize(read_count, 0ul);
    hist.push_back(frequency);
    prev_read_count = read_count;
  }
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
variance_bootstraps(const bool VERBOSE, const vector<double> &orig_hist,
		    const size_t bootstraps, const size_t orig_max_terms, 
		    const int diagonal, const double step_size, 
		    const double max_extrapolation, const double dupl_level,
		    const double tolerance, const size_t max_iter,
		    vector<double> &Ylevel_estimates,
		    vector< vector<double> > &yield_bootstrap_estimates) {
  // clear returning vectors
  yield_bootstrap_estimates.clear();
  
  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand()); 
  
  const double vals_size = accumulate(orig_hist.begin(), orig_hist.end(), 0.0);
  double vals_sum = 0.0;
  for(size_t i = 0; i < orig_hist.size(); i++)
    vals_sum += i*orig_hist[i];
  const double max_val = max_extrapolation/vals_sum;
  
  for (size_t iter = 0; 
       (iter < max_iter && yield_bootstrap_estimates.size() < bootstraps); ++iter) {
    
    vector<double> yield_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist, vals_sum, vals_size, hist);

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

    // interpolate complexity curve by random sampling w/out replacement
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
    // degree of approx is 1 less than max_terms
    max_terms = max_terms - (max_terms % 2 == 0);
    
    const ContinuedFractionApproximation 
      lower_cfa(diagonal, max_terms, step_size, max_extrapolation);
    //refit curve for lower bound
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
	yield_bootstrap_estimates.push_back(yield_vector);
	if (VERBOSE) cerr << '.';
	Ylevel_estimates.push_back(lower_cf.Ylevel(hist, dupl_level, vals_sum, 
						   max_val, tolerance, 100));
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
  if (yield_bootstrap_estimates.size() < bootstraps)
    throw SMITHLABException("too many iterations, poor sample");
}

static bool
single_estimates(const bool VERBOSE, vector<double> &hist,
		 size_t max_terms, const int diagonal, 
		 const double step_size, const double max_extrapolation, 
		 vector<double> &yield_estimate) {

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand()); 
  

  yield_estimate.clear();
  double vals_sum = 0.0;
  for(size_t i = 0; i < hist.size(); i++)
    vals_sum += i*hist[i];
  const double initial_distinct = accumulate(hist.begin(), hist.end(), 0.0);

  // const double max_val = max_extrapolation/vals_sum;
  // const double val_step = step_size/vals_sum;

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
    yield_estimate.push_back(sample_count_distinct(rng, umis, sample));
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
      yield_estimate.push_back(initial_distinct + t*lower_cf(t));
      sample_size += step_size;
    }
  }
  else{
    // FAIL!
    // lower_cf unacceptable, need to bootstrap to obtain estimates
    return false;
  }

  if (VERBOSE) {
    cerr << "CF_OFFSET_COEFF_ESTIMATES" << endl;
    copy(lower_cf.offset_coeffs.begin(), lower_cf.offset_coeffs.end(),
	 std::ostream_iterator<double>(cerr, "\n"));
    
    cerr << "CF_COEFF_ESTIMATES" << endl;
    copy(lower_cf.cf_coeffs.begin(), lower_cf.cf_coeffs.end(),
	 std::ostream_iterator<double>(cerr, "\n"));
  }

  // SUCCESS!!  
  return true;
}


static inline double
alpha_log_confint_multiplier(const double estimate,
			     const double variance, const double alpha) {
  const double inv_norm_alpha = gsl_cdf_ugaussian_Qinv(alpha/2.0);
  return exp(inv_norm_alpha*
	     sqrt(log(1.0 + variance/pow(estimate, 2))));
}


static void
vector_ci(const vector<vector<double> > &bootstrap_estimates,
	  const double ci_level, const vector<double> &yield_estimates,
	  vector<double> &lower_ci_lognormal, 
	  vector<double> &upper_ci_lognormal) {
  
  const double alpha = 1.0 - ci_level;
  assert(!bootstrap_estimates.empty());
  
  const size_t n_est = bootstrap_estimates.size();
  vector<double> estimates_row(bootstrap_estimates.size(), 0.0);
  for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {
    
    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];
    
    // sort to get confidence interval
    const double variance = gsl_stats_variance(&estimates_row[0], 1, n_est);
    const double confint_mltr = 
      alpha_log_confint_multiplier(yield_estimates[i], variance, alpha);
    lower_ci_lognormal.push_back(yield_estimates[i]/confint_mltr);
    upper_ci_lognormal.push_back(yield_estimates[i]*confint_mltr);
  }
}

static void
vector_median_ci(const vector<vector<double> > &bootstrap_estimates,
		 const double ci_level, vector<double> &yield_estimates,
		 vector<double> &lower_ci_lognormal, 
		 vector<double> &upper_ci_lognormal) {
  
  yield_estimates.clear();
  const double alpha = 1.0 - ci_level;
  assert(!bootstrap_estimates.empty());
  
  const size_t n_est = bootstrap_estimates.size();
  vector<double> estimates_row(bootstrap_estimates.size(), 0.0);
  for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {
    
    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];

    sort(estimates_row.begin(), estimates_row.end());
    const double median_estimate = 
      gsl_stats_median_from_sorted_data(&estimates_row[0], 1, n_est);
    
    // sort to get confidence interval
    const double variance = gsl_stats_variance(&estimates_row[0], 1, n_est);
    const double confint_mltr = 
      alpha_log_confint_multiplier(median_estimate, variance, alpha);
    yield_estimates.push_back(median_estimate);
    lower_ci_lognormal.push_back(median_estimate/confint_mltr);
    upper_ci_lognormal.push_back(median_estimate*confint_mltr);
  }
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
		      const vector<double> &yield_estimates,
		      const vector<double> &yield_lower_ci_lognormal,
		      const vector<double> &yield_upper_ci_lognormal) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOGNORMAL_LOWER_" << 100*c_level << "%CI\t"
      << "LOGNORMAL_UPPER_" << 100*c_level << "%CI" << endl;
  
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);
  
  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < yield_estimates.size(); ++i)
    out << (i + 1)*step_size << '\t' 
	<< yield_estimates[i] << '\t'
	<< yield_lower_ci_lognormal[i] << '\t' << yield_upper_ci_lognormal[i] << endl;
}



int
main(const int argc, const char **argv) {

  try {
    
    const size_t MIN_REQUIRED_COUNTS = 5;

    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 1000;
    double max_extrapolation = 1.0e10;
    double step_size = 1e6;
    size_t bootstraps = 100;
    int diagonal = -1;
    double c_level = 0.95;
    double tolerance = 1e-20;
    size_t max_iter = 100;
    double dupl_level = 0.5;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    bool SINGLE_ESTIMATE = false;
    
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
    opt_parse.add_opt("dupl_level",'d', "fraction of duplicate to predict "
		      "(default: " + toa(dupl_level) + ")",
		      false, dupl_level);
    opt_parse.add_opt("terms",'x',"maximum number of terms", false, 
		      orig_max_terms);
    //    opt_parse.add_opt("tol",'t', "numerical tolerance", false, tolerance);
    //    opt_parse.add_opt("max_iter",'i', "maximum number of iteration",
    //		      false, max_iter);
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
    opt_parse.add_opt("hist", 'H', 
		      "input is a text file containing the observed histogram",
		      false, HIST_INPUT);
    opt_parse.add_opt("quick",'Q',
		      "quick mode, estimate yield without bootstrapping for confidence intervals",
		      false, SINGLE_ESTIMATE);
    
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
    
    vector<double> counts_hist;
    double total_reads = 0.0;
    size_t max_observed_count = 0;
    double distinct_reads = 0.0;
    
    // LOAD VALUES
    if(HIST_INPUT){
      if(VERBOSE) 
	cerr << "INPUT_HIST" << endl;
      load_histogram(input_file_name, counts_hist);
      for(size_t i = 0; i < counts_hist.size(); i++){
	total_reads += i*counts_hist[i];
      }
      max_observed_count = counts_hist.size() - 1;
      distinct_reads = accumulate(counts_hist.begin(), counts_hist.end(), 0.0);
    }
    else{
      vector<double> values;
      if(VALS_INPUT){
	if(VERBOSE)
	  cerr << "VALS_INPUT" << endl;
	load_values(input_file_name, values);
      }
#ifdef HAVE_BAMTOOLS
      else if (BAM_FORMAT_INPUT && PAIRED_END){
	if(VERBOSE)
	  cerr << "PAIRED_END_BAM_INPUT" << endl;
	load_values_BAM_pe(input_file_name, values);
      }
      else if(BAM_FORMAT_INPUT){
	if(VERBOSE)
	  cerr << "BAM_INPUT" << endl;
	load_values_BAM_se(input_file_name, values);
      }
#endif
      else if(PAIRED_END){
	cerr << "PAIRED_END_BED_INPUT" << endl;
	load_values_BED_pe(input_file_name, values);
      }  
      else{
	if(VERBOSE)
	  cerr << "BED_INPUT" << endl;
	load_values_BED_se(input_file_name, values);
      }
    
    // JUST A SANITY CHECK
      total_reads = accumulate(values.begin(), values.end(), 0.0);


      max_observed_count = 
	static_cast<size_t>(*std::max_element(values.begin(), values.end()));

      distinct_reads = static_cast<double>(values.size());

    // BUILD THE HISTOGRAM
      counts_hist.clear();
      counts_hist.resize(max_observed_count + 1, 0.0);
      for (size_t i = 0; i < values.size(); ++i)
	++counts_hist[static_cast<size_t>(values[i])];
    }

    // for large initial experiments need to adjust step size
    // otherwise small relative steps do not account for variance
    // in extrapolation
    if(step_size < (total_reads/20)){
      step_size = std::max(step_size, step_size*round(total_reads/(20*step_size)));
      if(VERBOSE)
	cerr << "ADJUSTED_STEP_SIZE = " << step_size << endl;
    }
       
    
    const size_t distinct_counts = 
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
					bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
      cerr << "TOTAL READS     = " << total_reads << endl
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
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    }

    // catch if all reads are distinct or sample sufficiently deep
    if (max_observed_count < MIN_REQUIRED_COUNTS)
      throw SMITHLABException("sample not sufficiently deep or duplicates removed");
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // ESTIMATE COMPLEXITY CURVE

    if(VERBOSE)
      cerr << "[ESTIMATING YIELD CURVE]" << endl;
    vector<double> yield_estimates;
    bool SINGLE_ESTIMATE_SUCCESS = 
      single_estimates(VERBOSE, counts_hist, orig_max_terms, diagonal,
		       step_size, max_extrapolation, yield_estimates);

    if(SINGLE_ESTIMATE){
      // IF FAILURE, EXIT
      if(!SINGLE_ESTIMATE_SUCCESS)
	throw SMITHLABException("SINGLE ESTIMATE FAILED, NEED TO RUN FULL MODE FOR ESTIMATES");

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
  

      if(VERBOSE && !SINGLE_ESTIMATE_SUCCESS)
	cerr << "SINGLE ESTIMATE FAILED, NEED TO ESTIMATE MEDIAN FROM BOOTSTRAPS" << endl;
      if(max_iter == 0)
	max_iter = 4*bootstraps;
    
      vector<vector <double> > bootstrap_estimates;
      vector<double> Ylevel_estimates;
      variance_bootstraps(VERBOSE, counts_hist,  bootstraps, orig_max_terms,
			  diagonal, step_size, max_extrapolation, dupl_level,
			  tolerance, max_iter, Ylevel_estimates, bootstrap_estimates);
      
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
      if (VERBOSE)
	cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;

    // yield ci    
      vector<double> yield_upper_ci_lognormal, yield_lower_ci_lognormal,
	yield_upper_ci_quantile, yield_lower_ci_quantile;

      if(!SINGLE_ESTIMATE_SUCCESS){
	// use bootstrap estimates to obtain median estimates
	vector_median_ci(bootstrap_estimates, c_level, yield_estimates, 
			 yield_lower_ci_lognormal, yield_upper_ci_lognormal);
      }
      else{
	// use single estimates as the expected complexity curve
	vector_ci(bootstrap_estimates, c_level, yield_estimates, 
		  yield_lower_ci_lognormal, yield_upper_ci_lognormal);
      }
      
    // Y50 median and ci
      double Ylevel_median = 0.0;
      double Ylevel_lower_ci = 0.0;
      double Ylevel_upper_ci = 0.0;
      median_and_ci(Ylevel_estimates, c_level, Ylevel_median,
		    Ylevel_lower_ci, Ylevel_upper_ci);

    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
      if (VERBOSE) 
	cerr << "[WRITING OUTPUT]" << endl;
    
      write_predicted_curve(outfile, total_reads, c_level, step_size,
			    yield_estimates, yield_lower_ci_lognormal, 
			    yield_upper_ci_lognormal);

      if(VERBOSE){
	cerr << "Y" << 100*dupl_level << "MEASURE OF LIBRARY QUALITY: "
	   << "EXPECTED # READS TO REACH "
	     << 100*dupl_level << "% DUPLICATES" << endl;
	cerr << "Y" << 100*dupl_level << " = " << Ylevel_median << endl;
	cerr << 100*c_level << "%CI: (" << Ylevel_lower_ci << ", " 
	     << Ylevel_upper_ci << ")" << endl;
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
