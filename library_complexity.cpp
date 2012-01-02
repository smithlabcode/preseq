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


static void 
compute_num_and_denom_coeffs(const bool VERBOSE, const size_t max_terms, 
			     const vector<double> &initial_coeffs,
			     const double max_time, const double time_step,
			     const double tolerance, const double defect_tolerance,
			     vector<double> &contfrac_coeffs,
			     vector<double> &num_approx, vector<double> &denom_approx) {
  
  vector<double> coeffs(initial_coeffs);
  
  size_t max_number_approxs = static_cast<size_t>(round((max_terms-8)/2)); 
  
  // DETECT DEFECT, MOVE UP TABLE TO MORE CONSERVATIVE APPROX
  size_t denom_terms = max_terms/2;
  bool FOUND_APPROX = false;
  for (size_t i = 0; i < max_number_approxs && !FOUND_APPROX; i++) { 
    num_approx.clear(); 
    denom_approx.clear();
    bool DEFECT_FLAG = false;
    compute_pade_curve(VERBOSE, coeffs, max_time, time_step, 
		       defect_tolerance, tolerance, denom_terms, 
		       num_approx, denom_approx, DEFECT_FLAG);
    
    product_difference_alg_for_contfrac_coeffs(coeffs, 10, 
 					       contfrac_coeffs);
    quotient_difference_alg_for_contfrac_coeffs(coeffs, 20, 
						contfrac_coeffs);
    
    FOUND_APPROX = !DEFECT_FLAG;
    coeffs.pop_back();
    coeffs.pop_back();
    --denom_terms;
  }
  if (!FOUND_APPROX)
    throw SMITHLABException("no acceptable approximation found");
}


static void
compute_coverage(const bool VERBOSE, const vector<double> &counts_histogram,
                 const double max_time, const double time_step,
                 const double tolerance, const double defect_tolerance,
		 const size_t initial_max_terms, vector<double> &estimates) {
  
  // need max_terms (L+M+1) odd for convergence from below
  const size_t max_terms = initial_max_terms - (initial_max_terms % 2 == 1);
  
  size_t vals_sum = 0;
  for (size_t i = 0; i < counts_histogram.size(); i++)
    vals_sum += counts_histogram[i]*i;
  
  vector<double> coeffs(max_terms, 0.0);
  for (size_t j = 0; j < max_terms; j++)
    coeffs[j] = counts_histogram[j + 1]*pow(-1, j + 2)*(j + 1)/vals_sum;

  vector<double> num_approx, denom_approx; 
  vector<double> contfrac_coeffs;
  compute_num_and_denom_coeffs(VERBOSE, max_terms, coeffs, max_time, time_step,
			       tolerance, defect_tolerance, contfrac_coeffs,
			       num_approx, denom_approx);
  
  estimates.clear();
  for (size_t i = 0; i < num_approx.size(); ++i)
    estimates.push_back(std::min(1.0, 1 - num_approx[i]/denom_approx[i]));
}


static void
compute_distinct(const bool VERBOSE, const vector<double> &counts_histogram,
                 const double max_time, const double time_step,
                 const double tolerance, const double defect_tolerance,
		 const size_t initial_max_terms, vector<double> &estimates) {
  
  //need max_terms = L+M+1 to be even so that L+M is odd and we get convergence from above
  const size_t max_terms = initial_max_terms - (initial_max_terms % 2 == 0);
  
  vector<double> coeffs(max_terms, 0.0);
  for (size_t j = 0; j < max_terms; j++)
    coeffs[j] = counts_histogram[j + 1]*pow(-1, j+2);
  
  vector<double> num_approx, denom_approx; 
  vector<double> contfrac_coeffs;
  compute_num_and_denom_coeffs(VERBOSE, max_terms, coeffs, max_time, time_step,
			       tolerance, defect_tolerance, contfrac_coeffs,
			       num_approx, denom_approx);

  for (size_t i = 0; i < contfrac_coeffs.size(); ++i)
    cerr << std::setw(12) << std::fixed << std::setprecision(2) << contfrac_coeffs[i] << "\t" 
	 << std::setw(12) << std::fixed << std::setprecision(2) << coeffs[i] << endl;
  
  const double values_size = 
    accumulate(counts_histogram.begin(), counts_histogram.end(), 0.0);
  
  estimates.clear();
  double prev_estimate = 0.0;
  double t = 0.0;
  for (size_t i = 0; i < num_approx.size(); ++i, t += time_step) {
    const double curr_estimate = t*num_approx[i]/denom_approx[i] + values_size;
    estimates.push_back(std::max(prev_estimate, curr_estimate));
    prev_estimate = estimates.back();
  }
}


int
main(const int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    
    size_t max_terms = 1000;
    double tolerance = 1e-20;
    double max_time = 10;
    double time_step = 0.1;
    double defect_tolerance = 0.1;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool coverage = false;
    bool SMOOTH_HISTOGRAM = false; 
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(argv[0], "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("time",'m',"maximum time", false, max_time);
    opt_parse.add_opt("step",'s',"time step", false, time_step);
    opt_parse.add_opt("terms",'t',"maximum number of terms", false, max_terms);
    opt_parse.add_opt("tol", '\0', "general numerical tolerance",
		      false, tolerance);
    opt_parse.add_opt("def",'\0',"defect tolerance [for: abs(zero - pole)]",
		      false, defect_tolerance);
    opt_parse.add_opt("smooth",'\0',"smooth histogram (default: no smoothing)",
                      false, SMOOTH_HISTOGRAM);
    opt_parse.add_opt("coverage", '\0', "estimate genome coverage "
		      "(default: estimate distinct)", false, coverage);
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
    
    const size_t max_observed_count = 
      *std::max_element(values.begin(), values.end());
    
    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE: max_terms selection
    // needs to be dependent on max time heuristics suggest the method
    // is accurate even for large numbers, set
    // (max_time)^max_term<10^250.
    max_terms = std::min(max_terms, max_observed_count);
    if (max_terms % 2 == 0)
      --max_terms; 
    
    // BUILD THE HISTOGRAM
    vector<double> counts_histogram(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_histogram[values[i]];
    
    const size_t distinct_counts = 
      std::count_if(counts_histogram.begin(), counts_histogram.end(),
		    bind2nd(std::greater<size_t>(), 0));
    
    if (SMOOTH_HISTOGRAM) {
      const size_t smoothing_bandwidth = 4;
      const double smoothing_decay_factor = 15.0;
      smooth_histogram(smoothing_bandwidth, smoothing_decay_factor, counts_histogram);
    }
    
    if (VERBOSE)
      cerr << "TOTAL READS     = " << read_locations.size() << endl
	   << "DISTINCT READS  = " << distinct_reads << endl
	   << "DISTINCT COUNTS = " << distinct_counts << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_histogram[1] << endl
	   << "LIB SIZE GUESS  = " 
	   << std::fixed << std::setprecision(1) 
	   << distinct_reads/(1.0 - counts_histogram[1]/vals_sum) << endl;
    
    vector<double> estimates;
    if (coverage)
      compute_coverage(VERBOSE, counts_histogram, max_time, time_step, 
		       tolerance, defect_tolerance, max_terms, estimates);
    else compute_distinct(VERBOSE, counts_histogram, max_time, time_step, 
			  tolerance, defect_tolerance, max_terms, estimates);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    double time = 0.0;
    for (size_t i = 0; i < estimates.size(); ++i, time += time_step)
      out << std::fixed << std::setprecision(1) 
	  << (time + 1.0)*vals_sum << '\t' << estimates[i] << endl;
    
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
