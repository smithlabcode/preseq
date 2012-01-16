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
using std::ostream;

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
get_counts(const vector<SimpleGenomicRegion> &reads, vector<size_t> &counts) {
  counts.push_back(1);
  for (size_t i = 1; i < reads.size(); ++i)
    if (reads[i] == reads[i - 1]) counts.back()++;
    else counts.push_back(1);
}


// Checks if estimates are stable (derivative large) for the
// particular approximation (degrees of num and denom) at a specific
// point
bool
stable_estimate(const double t,
		const double prev_estim, 
		const double step_size,
		const double init_distinct,
		ContFracApprox &CFestimate){
  // stable if d/dt f(t) < vals_sum and f(t) <= prev_val+sample_per_time_step
  const double distinct_per_step = init_distinct*step_size;
  bool IS_STABLE = false;
  const double estimate = CFestimate.cont_frac_estimate.evaluate(t, CFestimate.get_depth());
  const double deriv = CFestimate.cont_frac_estimate.complex_deriv(t, CFestimate.get_depth());
  // we are using an CF approx that acts like x in limit
  // derivative must be positive and be less than the initial derivative
  if (deriv <= init_distinct && deriv >= 0.0)
    IS_STABLE = true;
  // make sure that the estimate is increasing in the time_step and is
  // below the initial distinct per step_size
  if (IS_STABLE && estimate >= prev_estim && estimate <= prev_estim + distinct_per_step)
    return IS_STABLE;
  else return false;
}

// Extrapolates the curve, for given values (step & max) and numbers
// of terms
static void
extrapolate_distinct(const bool VERBOSE, 
		     const vector<double> &counts_histogram,
		     const double max_value, const double step_size,
		     size_t max_terms, vector<double> &estimates) {
  
  // ensure that we will use an underestimate
  max_terms = max_terms - (max_terms % 2 == 1);
  
  // hist_sum = number of distinct in initial sample
  const double hist_sum = 
    accumulate(counts_histogram.begin(), counts_histogram.end(), 0.0);
  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_histogram.size(); i++)
    counts_sum += i*counts_histogram[i];
  
  vector<double> coeffs(max_terms, 0.0);
  for (size_t j = 0; j < max_terms; j++)
    coeffs[j] = counts_histogram[j + 1]*pow(-1, j+2);

  // initialize the cont frac estimator
  cont_frac cont_frac_estim(coeffs, 0, 0);
  ContFracApprox CFestimate(cont_frac_estim, max_terms);
  
  estimates.clear();
  while (max_terms > CFestimate.MINIMUM_ALLOWED_DEGREE && estimates.empty()) {
    estimates.push_back(hist_sum);
    double value = step_size;
    bool STABLE_ESTIMATE = true;
    while (value <= max_value && STABLE_ESTIMATE) {

      STABLE_ESTIMATE = STABLE_ESTIMATE &&
	stable_estimate(value, estimates.back(), 
			step_size, hist_sum, CFestimate);
      if (STABLE_ESTIMATE)
        estimates.push_back(hist_sum + CFestimate.cont_frac_estimate.evaluate(value, CFestimate.get_depth()));
      else {
	// estimates are unacceptable, move down in order
        estimates.clear();
        max_terms -= 2;
	CFestimate.set_depth(max_terms); 
      }
      value += step_size;
    }

    // if estimates.size() > 0 then the estimates are acceptable
    // output coeffs
    if (estimates.size() && VERBOSE) {
      vector<double> contfrac_coeffs, off_coeffs;
      off_coeffs = CFestimate.cont_frac_estimate.offset_coeffs;
      contfrac_coeffs = CFestimate.cont_frac_estimate.cf_coeffs;
      for(size_t i = 0; i < off_coeffs.size(); ++i)
	cerr << setw(12) << fixed << setprecision(2) 
	     << off_coeffs[i] << "\t" 
	     << setw(12) << fixed << setprecision(2) << coeffs[i] << endl;
      for (size_t i = 0; i < contfrac_coeffs.size(); ++i)
	cerr << setw(12) << fixed << setprecision(2) 
	     << contfrac_coeffs[i] << "\t" 
	     << setw(12) << fixed << setprecision(2) 
	     << coeffs[i+off_coeffs.size()] << endl;
    }
  }
}

int
main(const int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    string stats_outfile;

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
    
    // initialize ContFracApprox 
    if (SMOOTH_HISTOGRAM) 
      smooth_histogram(smoothing_bandwidth, 
		       smoothing_decay_factor, counts_histogram);
    
    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 0;
    while (counts_before_first_zero < counts_histogram.size() && 
	   counts_histogram[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
    max_terms = std::min(max_terms, counts_before_first_zero);     
    
    if (VERBOSE)
      cerr << "TOTAL READS     = " << read_locations.size() << endl
	   << "DISTINCT COUNTS = " << distinct_counts << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_histogram[1] << endl;
    
    vector<double> estimates;
    extrapolate_distinct(VERBOSE, counts_histogram, max_val,
			 val_step, max_terms, estimates);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    double val = 0.0;
    for (size_t i = 0; i < estimates.size(); ++i, val += val_step)
      out << std::fixed << std::setprecision(1) 
	  << (val + 1.0)*values_sum << '\t' << estimates[i] << endl;
    
    if (VERBOSE || !stats_outfile.empty()) {
      // TODO: WITH BETTER OUTPUT
      std::ofstream stats_of;
      if (!stats_outfile.empty()) stats_of.open(stats_outfile.c_str());
      ostream stats_out(stats_outfile.empty() ? cerr.rdbuf() : stats_of.rdbuf());
      
      const double upper_bound = 
	upperbound_librarysize(counts_histogram, max_terms)+distinct_counts;
      stats_out << "Chao87_lower_bound" << "\t" 
		<< chao87_lowerbound_librarysize(counts_histogram) << endl;
      stats_out << "ChaoLee92_lower_bound" << "\t" 
		<< cl92_lowerbound_librarysize(counts_histogram) << endl;
      stats_out << "cont_frac_lower_bound" << "\t"
		<< lowerbound_librarysize(counts_histogram, upper_bound, 
					  val_step, max_val, max_terms) << endl;
      stats_out << "cont_frac_upper_bound" << "\t" << upper_bound << endl;
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
