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

#include <fstream>
#include <numeric>
#include <vector>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;

using std::setw;
using std::fixed;
using std::setprecision;

using smithlab::log_sum_log_vec;


static inline double
weight_exponential(const double dist, double decay_factor) {
  return std::pow(0.5, decay_factor*dist);
}


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
  updated_hist.swap(hist);
}


static void
get_counts(const vector<SimpleGenomicRegion> &reads,
           vector<size_t> &counts) {
  counts.push_back(1);
  for (size_t i = 1; i < reads.size(); ++i)
    if (reads[i] == reads[i - 1]) counts.back()++;
    else counts.push_back(1);
}


// Checks if estimates are stable (derivative large) for the
// particular approximation (degrees of num and denom) at a specific
// point
bool
cont_frac_distinct_stable(const bool VERBOSE, const double val,
                          const double prev_estim, 
			  const double step_size,
			  const double init_distinct,
                          ContFracApprox &CFestimate){
  // stable if d/dt f(t) < vals_sum and f(t) <= prev_val+sample_per_time_step
  const double distinct_per_step = init_disinct*step_size;
  bool IS_STABLE = false;
  const double estimate = CFestimate.cont_frac_estimate.evaluate(val);
  const double deriv =  
    = CFestimate.cont_frac_estimate.complex_deriv(val);
  // we are using an CF approx that acts like x in limit
  // derivative must be positive and be less than the initial derivative
  if(deriv <= init_distinct && deriv >= 0.0){
    IS_STABLE = true;
  }
  //make sure that the estimate is increasing in the time_step
  // and is below the initial distinct per step_size
  if (IS_STABLE && estimate >= prev_estim &&
      current_val <= prev_estim + distinct_per_step) {
    return IS_STABLE;  //estimate is stable, exit_success
  }
  else {
    IS_STABLE = false;
    if (VERBOSE) {
      cerr << "error_loc" << "\t" << "estimate" << "\t" 
	   << "deriv" << endl;
      cerr << val << "\t" <<  CFestimate.cont_frac_estimate.evaluate(val)
	   << "\t" << CFestimate.cont_frac_estimate.complex_deriv(val)
	   << endl;
    } 
// if estimate is unstable, exit_failure, else exit_success
    return IS_STABLE; 
  }
}

// Extrapolates the curve, for given values (step & max) and numbers
// of terms
static void
extrapolate_distinct(const bool VERBOSE, 
		     const vector<double> &counts_histogram,
		     const double max_value, const double step_size,
		     const double tolerance, const double deriv_delta, 
		     const size_t initial_max_terms, 
		     vector<double> &estimates) {
  // ensure that we will use an underestimate
  size_t max_terms = initial_max_terms - (initial_max_terms % 2 == 1);
  
  // hist_sum = number of distinct in initial sample
  const double hist_sum = 
  accumulate(counts_histogram.begin(), counts_histogram.end(), 0.0);
  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_histogram.size(); i++)
    vals_sum += i*counts_histogram[i];
  
  vector<double> coeffs(max_terms, 0.0);
  for (size_t j = 0; j < max_terms; j++)
    coeffs[j] = counts_histogram[j + 1]*pow(-1, j+2);

  cont_frac cont_frac_estim(coeffs, 0, 0);
  ContFracApprox CFestimate(cont_frac_estim);
  
  estimates.clear();
  while (max_terms > CFestimate.MINIMUM_ALLOWED_DEGREE &&
	 estimates.size() == 0) {
    estimates.push_back(hist_sum);
    double value = step_size;
    bool STABLE_ESTIMATE = true;
    while (value <= max_value && STABLE_ESTIMATE) {

      STABLE_ESTIMATE = STABLE_ESTIMATE &&
	cont_frac_distinct_stable(VERBOSE, val, prev_estim, 
				  step_size, hist_sum, CFestimate);
      if(STABLE_ESTIMATE)
        estimates.push_back(values_size 
			    + CFestimate.cont_frac_estimate.cf_approx(val));
      else{
        estimates.clear();
        max_terms -= 2;
      }
      val += step_size;
    }
    if (estimates.size()) { //if estimates.size() > 0
      if (VERBOSE) {
	vector<double> contfrac_coeffs;
        contfrac_estimate.get_cf_coeffs(contfrac_coeffs);
        for (size_t i = 0; i < contfrac_coeffs.size(); ++i)
          cerr << setw(12) << fixed << setprecision(2) 
	       << contfrac_coeffs[i] << "\t" 
	       << setw(12) << fixed << setprecision(2) << coeffs[i] << endl;
      }
    }
    else{
      coeffs.push_back();
      coeffs.push_back();
      ContFracApprox.set_ps_coeffs(coeffs);
    }
  }
}


static double
chao87_lowerbound_librarysize(const vector<double> &counts_histogram){
  return(accumulate(counts_histogram.begin(), counts_histogram.end(), 0.0)
         + counts_histogram[1]*counts_histogram[1]/(2*counts_histogram[2]));
}


//Chao & Lee (JASA 1992) lower bound
static double 
chao_lee_librarysize(const vector<double> &counts_histogram){ 
  double sample_size = 0.0;
  for(size_t i = 0; i <  counts_histogram.size(); i++)
    sample_size += counts_histogram[i]*i;
  const double distinct_vals = accumulate(counts_histogram.begin(), counts_histogram.end(), 0.0);
  const double estim_coverage = 1.0 - counts_histogram[1]/sample_size;
  const double naive_lowerbound = distinct_vals/estim_coverage;
  vector<double> log_coeff_var;
  for(size_t i = 2; i < counts_histogram.size(); i++){
    if(counts_histogram[i] > 0)
      log_coeff_var.push_back(log(i) + log(i - 1) + log(counts_histogram[i]));
  }
  const double coeff_var = 
    naive_lowerbound*exp(log_sum_log_vec(log_coeff_var, log_coeff_var.size())
			 -log(sample_size) - log(sample_size - 1)) - 1;
  return(naive_lowerbound + counts_histogram[1]*coeff_var/estim_coverage);
}


//// TODO: good comment here!!!!!!!!!!!!!
static double
upperbound_librarysize(const bool VERBOSE, const vector<double> &counts_histogram,
                       const size_t initial_max_terms){

  //need max_terms = L+M+1 to be even so that L+M is odd so that we can take lim_{t \to \infty} [L+1, M]
  size_t max_terms = initial_max_terms - (initial_max_terms % 2 == 1); 
  vector<double> coeffs(max_terms, 0.0);
  for(size_t j = 0; j < max_terms; j++)
    coeffs[j] = counts_histogram[j + 1]*pow(-1, j+2);
  
  vector<double> denom_vec;
  vector<double> num_vec;
  // TOD: MAGIC!!!!!!!!!!!!!!!
  while(max_terms >= 12){
    const size_t numer_size = max_terms/2;  //numer_size = L+1, denom_size = M
    const size_t denom_size = max_terms - numer_size;
    if(numer_size != denom_size)
      cerr << "num size = " << numer_size << ", denom size = " << denom_size << "\n";
    assert(numer_size == denom_size);
    bool ACCEPT_APPROX = compute_pade_coeffs(coeffs, numer_size, denom_size, num_vec, denom_vec); 
    if(ACCEPT_APPROX && num_vec.back()/denom_vec.back() > 0){
      if(VERBOSE){
        cerr << "numerator coeffs = ";
        for(size_t i = 0; i < num_vec.size(); i++)
          cerr << num_vec[i] << ", ";
        cerr << "\ndenomimator coeffs = ";
        for(size_t i = 0; i < denom_vec.size(); i++)
          cerr << denom_vec[i] << ", ";
        cerr << "\n";
      }
      break;
    }
    else{
      if(VERBOSE)
        cerr << "unacceptable approx, number of terms = " << max_terms << "\n";
      denom_vec.clear();
      num_vec.clear();
      max_terms -= 2;
    }
  }
  return(num_vec.back()/denom_vec.back());
}

static double
lowerbound_librarysize(const bool VERBOSE, const vector<double> &counts_histogram,
                        const double upper_bound, //use pade_lowerbound_librarysize(.) to compute
                       const double time_step, const double max_time,
                       const double deriv_delta, const double tolerance,
                       const size_t initial_max_terms){ 
  double vals_sum = 0.0;
  for(size_t i = 0; i < counts_histogram.size(); i++)
    vals_sum += i*counts_histogram[i];
  
  const double distinct_vals = accumulate(counts_histogram.begin(), counts_histogram.end(), 0.0); 
  
  size_t max_terms = initial_max_terms - (initial_max_terms % 2 == 0);
  
  vector<double> contfrac_coeffs;
  vector<double> contfrac_offsetcoeffs;
  cont_frac contfrac_estimate(contfrac_coeffs, contfrac_offsetcoeffs, 2, 0);
  
  vector<double> coeffs(max_terms, 0.0);
  for(size_t j = 0; j < max_terms; j++)
    coeffs[j] = counts_histogram[j+1]*pow(-1, j+2);
  
  vector<double> possible_maxima_loc;
  vector<double> possible_maxima;
  // MAGIC!!!!!!!!!!!!!!!!!!!
  while (max_terms > 6){
    contfrac_estimate.compute_cf_coeffs(coeffs, max_terms);
    // TODO: "TIME"
    double t = time_step;
    double prev_deriv = 0.0;
    double current_deriv = contfrac_estimate.cf_deriv_complex(t, deriv_delta, tolerance);
    double current_val = distinct_vals;
    double prev_val = distinct_vals;
    while(t < max_time){
      current_deriv = contfrac_estimate.cf_deriv_complex(t, deriv_delta, tolerance);
      current_val = contfrac_estimate.cf_approx(t, tolerance)+distinct_vals;
      //if derivative or estimate is not acceptable, choose different order approx
      if(fabs(current_deriv) > vals_sum || current_val > upper_bound){
	//so that we can know we need to go down in approx
        possible_maxima_loc.clear();  
        possible_maxima.clear();
	// TODO: NO BREAK!!!!!!!!!!!!!!
        break; //out of t loop
      }
      if(current_deriv*prev_deriv < 0.0 && current_val < upper_bound){
        possible_maxima_loc.push_back(contfrac_estimate.locate_zero_cf_deriv(t, t-time_step,
                                                                             deriv_delta, tolerance));
        possible_maxima.push_back(contfrac_estimate.cf_approx(possible_maxima_loc.back(), tolerance)+distinct_vals);
      }
      else if (current_val < prev_val && prev_val < upper_bound) {
        possible_maxima_loc.push_back(t-time_step);
        possible_maxima.push_back(prev_val);
      }
      prev_deriv = current_deriv;
      prev_val = current_val;
      t += time_step;
    }
    if (possible_maxima.size() > 0) {
      break; //out of max_terms loop
    }
    else{
      const vector<double> void_contfrac_coeffs;
      const vector<double> void_offset_coeffs;
      contfrac_estimate.set_cf_coeffs(void_contfrac_coeffs);
      contfrac_estimate.set_offset_coeffs(void_offset_coeffs);
      max_terms -= 2;
    }
  }
  
  possible_maxima_loc.push_back(max_time); //include boundary
  possible_maxima.push_back(contfrac_estimate.cf_approx(max_time, tolerance)+distinct_vals);

  /// TOD: BETTER VERBOSE COMMENTS (PARSE-ABLE)
  if (VERBOSE) {
    cerr << "possible maxima = ";
    for(size_t i = 0; i < possible_maxima_loc.size(); i++)
      cerr << possible_maxima_loc[i] << ", ";
    cerr << "\nvalues = ";
    for(size_t i = 0; i < possible_maxima.size(); i++)
      cerr << possible_maxima[i] << ", ";
    cerr << "\n";
    cerr << "upper bound = " << upper_bound << "\n";
  }
  double current_global_max = 0.0;
  double current_global_max_loc = 0.0;
  for(size_t j = 0; j < possible_maxima_loc.size(); j++){
    const double test_val = possible_maxima[j];
    if(test_val > current_global_max && test_val < upper_bound ){ 
      current_global_max = test_val;
      current_global_max_loc = possible_maxima_loc[j];
    }
  }
  if(VERBOSE)
    cerr << "CHOSEN_GLOBAL_MAX=" << current_global_max << "\t"
	 << "LOC=" << current_global_max_loc << endl;
  return current_global_max;
}

/* Chao & Lee (JASA 1992) lower bound
 */
static double 
cl92_library_size_lowerbound(const vector<double> &vals_hist) {

  double sample_size = 0.0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    sample_size += i*vals_hist[i]; 
  const double distinct = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
  
  // The "coverage" below is an estimate??
  const double coverage = 1.0 - vals_hist[1]/sample_size;
  
  const double naive_lowerbound = distinct/coverage;

  /// This is sum(i(i-1)f_i) from 2.12 and 2.13 of CL92 paper
  vector<double> log_cv_terms;
  for (size_t i = 2; i < vals_hist.size(); ++i) 
    if (vals_hist[i] > 0)
      log_cv_terms.push_back((log(i) + log(i - 1)) + log(vals_hist[i]));
  const double cv_terms_sum = 
    (log_cv_terms.empty()) ? 0.0 : 
    exp(log_sum_log_vec(log_cv_terms, log_cv_terms.size()));
  
  const double coeff_variation_est = 
    max(naive_lowerbound*cv_terms_sum/
	(sample_size*(sample_size - 1.0)) - 1.0, 0.0);
  
  // corrected coefficient of variation below is from 2.13 in
  // Chao & Lee (JASA 1992)
  const double cv_bias_correction = 1.0 + 
    ((1.0 - coverage)*cv_terms_sum)/((sample_size - 1.0)*coverage);
  
  const double coeff_variation_corr = 
    max(coeff_variation_est*cv_bias_correction,	0.0);
  
  return naive_lowerbound + 
    coeff_variation_corr*(sample_size*(1.0 - coverage)/coverage);
}


//Chao & Lee (JASA 1992) lower bound
double 
chao_lee_lowerbound_librarysize(const vector<double> &vals_hist){ 
  double sample_size = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++)
    sample_size += i*vals_hist[i]; 
  const double distinct = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
  double coverage = 1 - vals_hist[1]/sample_size;
  double naive_lowerbound = distinct/coverage;
  vector<double> log_cv_terms;
  for(size_t i = 2; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0) {
      log_cv_terms.push_back(log(naive_lowerbound) + log(i) + log(i-1)
                             -log(sample_size) - log(sample_size-1));
    }
  }
  double coeff_variation = 0.0;
  if(log_cv_terms.size() > 0)
    coeff_variation = max(exp(log_sum_log_vec(log_cv_terms, log_cv_terms.size()))-1, 0.0);
  
  for(size_t i = 0; i < log_cv_terms.size(); i++)
    log_cv_terms[i] -= log(naive_lowerbound);
  const double corrected_coeff_variation = 
    coeff_variation*(1+sample_size*(1-coverage)
		     *exp(log_sum_log_vec(log_cv_terms, 
					  log_cv_terms.size()))/coverage);
  
  return naive_lowerbound + sample_size*(1-coverage)*corrected_coeff_variation/coverage;
}


int
main(const int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    stirng stats_outfile;

    size_t max_terms = 1000;
    double tolerance = 1e-20;
    double max_extrapolation = 1e10;
    double step_size = 1e6;
    double deriv_delta = 1e-8;
    size_t smoothing_bandwidth = 4;
    double smoothing_decay_factor = 15.0;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool SMOOTH_HISTOGRAM = false; 
    bool LIBRARY_SIZE = false;
    
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
    opt_parse.add_opt("terms",'t',"maximum number of terms", false, max_terms);
    opt_parse.add_opt("tol", '\0', "general numerical tolerance",
		      false, tolerance);
    opt_parse.add_opt("delta", '\0', "derivative step size",
                      false, deriv_delta);
    opt_parse.add_opt("smooth",'\0',"smooth histogram (default: no smoothing)",
                      false, SMOOTH_HISTOGRAM);
    opt_parse.add_opt("bandwidth", '\0', "smoothing bandwidth",
                      false, smoothing_bandwidth);
    opt_parse.add_opt("decay", '\0', "smoothing decay factor",
                      false, smoothing_decay_factor);
    opt_parse.add_opt("library_size", '\0', "estimate library size "
		      "(default: estimate distinct)", false, LIBRARY_SIZE);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
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
    ReadBEDFile(input_file_name, read_locations);
    if (!check_sorted(read_locations))
      throw SMITHLABException("read_locations not sorted");
    
    // OBTAIN THE COUNTS FOR DISTINCT READS
    vector<size_t> values;
    get_counts(read_locations, values);
    
    const size_t distinct_reads = values.size();
    
    // JUST A SANITY CHECK
    const size_t n_reads = read_locations.size();
    const size_t vals_sum = accumulate(values.begin(), values.end(), 0ul);
    assert(vals_sum == n_reads);
    
    const double max_time = max_extrapolation/static_cast<double>(vals_sum);
    const double time_step = step_size/static_cast<double>(vals_sum);
    
    const size_t max_observed_count = 
      *std::max_element(values.begin(), values.end());

        
    // BUILD THE HISTOGRAM
    vector<double> counts_histogram(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_histogram[values[i]];
    
    const size_t distinct_counts = 
      std::count_if(counts_histogram.begin(), counts_histogram.end(),
		    bind2nd(std::greater<size_t>(), 0));
          
    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t max_count_b4_zero = 0;
    size_t indx = 1;
    while(counts_histogram[indx] > 0){
      max_count_b4_zero++;
      indx++;
    }
    max_terms = std::min(max_terms, max_count_b4_zero);     
    
    if (SMOOTH_HISTOGRAM) 
      smooth_histogram(smoothing_bandwidth, smoothing_decay_factor, counts_histogram);
    
    
    if (VERBOSE)
      cerr << "TOTAL READS     = " << read_locations.size() << endl
	   << "DISTINCT READS  = " << distinct_reads << endl
	   << "DISTINCT COUNTS = " << distinct_counts << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_histogram[1] << endl;
    
    vector<double> estimates;
    if (LIBRARY_SIZE){
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      
      out << "Chao 1987 lower bound" << "\t" 
	  << chao87_lowerbound_librarysize(counts_histogram) << endl;
      out << "Chao-Lee estimate" << "\t" 
	  << chao_lee_librarysize(counts_histogram) << endl;
      const double upper_bound = upperbound_librarysize(VERBOSE, counts_histogram, max_terms)+distinct_reads;
      out << "Continued Fraction lower bound" << "\t"
	  << lowerbound_librarysize(VERBOSE, counts_histogram, upper_bound, time_step, max_time, 
				    deriv_delta, tolerance, max_terms) << endl;
      out << "Continued Fraction upper bound" << "\t" << upper_bound << endl;
    }
    else { 
      extrapolate_distinct(VERBOSE, counts_histogram, max_time, time_step,
			   tolerance, deriv_delta, max_terms, estimates);
      
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      
      double t = 0.0;
      for (size_t i = 0; i < estimates.size(); ++i, t += time_step)
        out << std::fixed << std::setprecision(1) 
	    << (t + 1.0)*vals_sum << '\t' << estimates[i] << endl;
    }

    if (VERBOSE) {
      // TODO: WITH BETTER OUTPUT
      cerr << chao_lee_lowerbound_librarysize(counts_histogram) << endl
	   << cl92_library_size_lowerbound(counts_histogram) << endl;
    }
    if (!stats_outfile.empty()) {
      // TOD: no recompute!!!!!!!
      std::ostream stats_out(stats_outfile.c_str());
      stats_out << chao_lee_lowerbound_librarysize(counts_histogram) << endl
		<< cl92_library_size_lowerbound(counts_histogram) << endl;
      
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
