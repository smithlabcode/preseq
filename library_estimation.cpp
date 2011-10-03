/*    library_estimation:
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *                       Timothy Daley
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


#include "my_pade.hpp"

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>


#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

using std::max;
using std::setw;
using std::fabs;
using std::ceil;
using std::greater;
using std::numeric_limits;


static inline bool
is_greater(const double val1, const double val2){
  return((val1 > val2));
}

static void
get_counts(const vector<SimpleGenomicRegion> &read_locations,
           vector<size_t> &values) {
  size_t count = 1;
  for (size_t i = 1; i < read_locations.size(); i++)
    if (read_locations[i] == read_locations[i-1])
      ++count;
    else {
      values.push_back(count);
      count = 1;
    }
}

void
compute_coverage(const vector<double> &vals_hist,
                 const double max_time,
                 const double time_step,
                 const double tolerance,
                 const double defect_error,
                 const bool VERBOSE,
                 vector<double> &estimate_vec,
                 size_t &max_terms){
  if(max_terms % 2 == 1)
    max_terms--; //need max_terms = L+M+1 to be odd for convergence from below
  
  size_t vals_sum = 0;
  for(size_t i = 0; i < vals_hist.size(); i++)
    vals_sum += vals_hist[i]*i;
  
  vector<double> numerator_approx;
  vector<double> denominator_approx; 
  vector<double> summands;
  vector<double> coeffs(max_terms, 0.0);
  for(size_t j = 0; j < max_terms; j++)
    coeffs[j] = vals_hist[j+1]*pow(-1, j+2)*(j+1)/vals_sum;
  
  size_t denom_size = max_terms/2; //denom_size = M
  
  bool defect_flag = false;
  size_t max_number_approxs = static_cast<size_t>(round((max_terms-8)/2)); //don't want to use too few terms
  for(size_t i = 0; i < max_number_approxs; i++){ //detect defect, move up table to more conservative approx
    compute_pade_curve(coeffs, max_time, time_step, defect_error, tolerance,
                       denom_size, VERBOSE, numerator_approx, denominator_approx, defect_flag);
    if(defect_flag == false)
      break;  //no allowable defect detected, accept approx
    coeffs.pop_back();
    coeffs.pop_back();
    denom_size--;
    numerator_approx.clear();
    denominator_approx.clear();
  }
  if(defect_flag == true){
    cerr << "no alowable approximation, decrease max_time to increase number of terms used \n";
  }
  else{
    estimate_vec.clear();
    for(size_t i = 0; i < numerator_approx.size(); i++)
      estimate_vec.push_back(std::min(1.0, 1 - numerator_approx[i]/denominator_approx[i]));
  }

}

void
compute_distinct(const vector<double> &vals_hist,
                 const double max_time,
                 const double time_step,
                 const double tolerance,
                 const double defect_error,
                 const bool VERBOSE,
                 vector<double> &estimate_vec,
                 size_t &max_terms){
  if(max_terms % 2 == 0)
    max_terms--; //need max_terms = L+M+1 to be even so that L+M is odd and we get convergence from above
  vector<double> summands;
  vector<double> coeffs(max_terms, 0.0);
  for(size_t j = 0; j < max_terms; j++)
    coeffs[j] = vals_hist[j+1]*pow(-1, j+2);
  
  const double values_size = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
  
  size_t denom_size = max_terms/2; //denom_size = M
  
  vector<double> numerator_approx;
  vector<double> denominator_approx;
  
  bool defect_flag = false;
  size_t max_number_approxs = static_cast<size_t>(round((max_terms-8)/2));
  for(size_t i = 0; i < max_number_approxs; i++){ //detect defect, move up table to more conservative approx
    compute_pade_curve(coeffs, max_time, time_step, defect_error, tolerance,
                       denom_size, VERBOSE, numerator_approx, denominator_approx, defect_flag);
    if(defect_flag == false)
      break;  //no allowable defect detected, accept approx
    coeffs.pop_back();
    coeffs.pop_back();
    denom_size--;
    numerator_approx.clear();
    denominator_approx.clear();
  }
  if(defect_flag == true){
    cerr << "no alowable approximation, decrease max_time to increase number of terms used \n";
  }
  else{
    double last_estimate = 0.0;
    estimate_vec.clear();
    double t = 0.0;
    for(size_t i = 0; i < numerator_approx.size(); i++){
      estimate_vec.push_back(std::max(last_estimate, t*numerator_approx[i]/denominator_approx[i] + values_size));
      t += time_step;
    }
    last_estimate = estimate_vec.back();
  }
}

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


int
main(int argc, const char **argv){


  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    double tolerance = 1e-20;
    size_t max_terms = 1000;
    double max_time = 10;
    double time_step = 0.1;
    double defect_error = time_step/10;
    bool coverage = false;
    bool genome_coverage = false;
    bool smooth_hist = false; 

    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("library_estimation", "",
                           "*.bed");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("time_step",'p',"time step",
                      false, time_step);
    opt_parse.add_opt("tolerance", 't', "Numerical tolerance",
                     false, tolerance);
    opt_parse.add_opt("max_time",'s',"maximum time",
                     false, max_time);
    opt_parse.add_opt("defect_error",'e',"allowable error between zero and pole in approximation",
                      false, defect_error);
    opt_parse.add_opt("smooth_hist",'h',"include to smooth histogram",
                      false, smooth_hist);
    opt_parse.add_opt("coverage",'\0', "estimate coverage, default estimates distinct",
                      false, coverage);
    opt_parse.add_opt("genome_coverage",'\0',"estimate genome coverage, default estimates distinct",
                      false, genome_coverage);
    




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
    
    vector<size_t> values;
    vector<SimpleGenomicRegion> read_locations;
    ReadBEDFile(input_file_name, read_locations);
    if (!check_sorted(read_locations))
      throw RMAPException("read_locations not sorted");
    
    get_counts(read_locations, values);
    
    const size_t values_size = values.size();
    const size_t vals_sum = accumulate(values.begin(), values.end(), 0);

    const size_t max_value = *std::max_element(values.begin(),values.end());
    vector<double> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i){
      ++vals_hist[static_cast<size_t>(values[i])];
    }
    
    if(smooth_hist){
      size_t smoothing_bandwidth = 4;
      double smoothing_decay_factor = 15.0;
      smooth_histogram(smoothing_bandwidth, smoothing_decay_factor, vals_hist);
    }
    
    vector<double> estimate_vec;
    
    //max_terms selection needs to be dependent on max time
    // heuristics suggest the method is accurate even for large numbers,
    // set (max_time)^max_term < 10^250
    max_terms = std::min(max_value, 
                         static_cast<size_t>(round(250*log(10)/log(max_time))));
    if(VERBOSE)
      cerr << "max term = " << max_terms << "\n";
    
    
    if(coverage){
      compute_coverage(vals_hist, max_time, time_step, tolerance,
                       defect_error, VERBOSE, estimate_vec, max_terms);
    }
    else if(genome_coverage){
      //do nothing
    }
    else{
      compute_distinct(vals_hist, max_time, time_step, tolerance,
                       defect_error, VERBOSE, estimate_vec, max_terms);
    }
    
    ostream* out = (outfile.empty()) ? 
    &std::cout : new std::ofstream(outfile.c_str());
    if(genome_coverage){
      *out << static_cast<double>(values_size)/(1-static_cast<double>(vals_hist[1])/static_cast<double>(vals_sum)) 
       << endl;
    }
    else{
      double time = 0.0;
      for(size_t i = 0; i < estimate_vec.size(); i++){
        *out << (time+1.0)*vals_sum << "\t" << estimate_vec[i] << endl;
        time += time_step;
      }
    }
    

  }  
  catch (RMAPException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
