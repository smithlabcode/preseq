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
#include <complex>

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

using std::complex;
using std::norm;
using std::real;
using std::imag;



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
  
/*
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
 */

bool
test_stability_cf_distinct(cont_frac cf_estimate,
                           const double time,
                           const double prev_val,
                           const double dx,
                           const double tolerance,
                           const double vals_sum,
                           const double samples_per_time_step,
                           const bool VERBOSE){
  // stable if d/dt f(time) < 1 and f(time) <= prev_val+sample_per_time_step
  bool return_val = false;
  const double current_val = cf_estimate.cf_approx(time, tolerance);
  double test_val = cf_estimate.cf_deriv_complex(time, dx, tolerance);
  if(test_val <= vals_sum && test_val >= 0.0){
    return_val = true;
  }
  if(return_val == true &&
     current_val <= prev_val + samples_per_time_step &&
     current_val >= prev_val){
    return(true);
  }
  else{
    if(VERBOSE){
    cerr << "error found at " << time << ", cf approx = " 
      << cf_estimate.cf_approx(time, tolerance) << ", deriv = "
      << cf_estimate.cf_deriv_complex(time, dx, tolerance) << 
      ", " << (cf_estimate.cf_approx(time+dx, tolerance)-cf_estimate.cf_approx(time, tolerance))/dx 
      << ", vals_sum = " << vals_sum << "\n";
    } 
    return(false);
  }
}
  
  
void
compute_distinct(const vector<double> &vals_hist,
                 const double vals_sum,
                 const double max_time,
                 const double time_step,
                 const double dx,
                 const double tolerance,
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
  
  while(max_terms > 10){
    vector<double> cf_coeffs;
    vector<double> offset_coeffs;
    cont_frac cf_estimate(cf_coeffs, offset_coeffs, 0, 0);
    cf_estimate.compute_cf_coeffs(coeffs, max_terms);
    estimate_vec.push_back(values_size);
    double time = time_step;
    while(time <= max_time){
      if(test_stability_cf_distinct(cf_estimate, time, estimate_vec.back()-values_size, dx,
                                    tolerance, vals_sum, vals_sum*time_step, VERBOSE)){
        estimate_vec.push_back(values_size + cf_estimate.cf_approx(time, tolerance));
      }
      else{
        if(VERBOSE)
          cerr << ", number of terms = " << max_terms << "\n";
        estimate_vec.clear();
        max_terms -= 2;
        break;
      }
      time += time_step;
    }
    if(estimate_vec.size() > 1){ //break from max_terms loop
      break;
    }
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

double 
chao_lee_lowerbound_librarysize(const vector<double> &vals_hist){ //Chao & Lee (JASA 1992) lower bound
  double sample_size = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++)
    sample_size += i*vals_hist[i]; 
  const double distinct = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
  double coverage = 1 - vals_hist[1]/sample_size;
  double naive_lowerbound = distinct/coverage;
  vector<double> log_cv_terms;
  for(size_t i = 2; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      log_cv_terms.push_back(log(naive_lowerbound) + log(i) + log(i-1)
                             -log(sample_size) - log(sample_size-1));
    }
  }
  double coeff_variation = 0.0;
  if(log_cv_terms.size() > 0)
    coeff_variation = max(exp(log_sum_log_vec(log_cv_terms, log_cv_terms.size()))-1, 0.0);
  
  for(size_t i = 0; i < log_cv_terms.size(); i++)
    log_cv_terms[i] -= log(naive_lowerbound);
  double corrected_coeff_variation = coeff_variation*(1+sample_size*(1-coverage)
                                                      *exp(log_sum_log_vec(log_cv_terms, 
                                                                           log_cv_terms.size()))/coverage);
  
  return(naive_lowerbound + sample_size*(1-coverage)*corrected_coeff_variation/coverage);
}

double
chao87_lowerbound_librarysize(const vector<double> &vals_hist){
  return(accumulate(vals_hist.begin(), vals_hist.end(), 0.0)
         + vals_hist[1]*vals_hist[1]/(2*vals_hist[2]));
}


double
pade_upperbound_librarysize(const vector<double> &vals_hist,
                            const bool VERBOSE,
                            size_t max_terms){

  vector<double> coeffs(max_terms, 0.0);
  for(size_t j = 0; j < max_terms; j++)
    coeffs[j] = vals_hist[j+1]*pow(-1, j+2);
    
  if(max_terms % 2 == 1)
    max_terms--; //need max_terms = L+M+1 to be even so that L+M is odd so that we can take lim_{t \to \infty} [L+1, M]
  vector<double> denom_vec;
  vector<double> num_vec;
  while(max_terms >= 12){
    size_t numer_size = max_terms/2;  //numer_size = L+1, denom_size = M
    size_t denom_size = max_terms - numer_size;
    if(numer_size != denom_size)
      cerr << "num size = " << numer_size << ", denom size = " << denom_size << "\n";
    assert(numer_size == denom_size);
    bool accept_approx = compute_pade_coeffs(coeffs, numer_size, denom_size, num_vec, denom_vec); 
    if(VERBOSE){
      cerr << "numerator coeffs = ";
      for(size_t i = 0; i < num_vec.size(); i++)
        cerr << num_vec[i] << ", ";
      cerr << "\ndenomimator coeffs = ";
      for(size_t i = 0; i < denom_vec.size(); i++)
        cerr << denom_vec[i] << ", ";
      cerr << "\n";
    }
    if(accept_approx && num_vec.back()/denom_vec.back() > 0)
      break;
    else{
      denom_vec.clear();
      num_vec.clear();
      max_terms -= 2;
    }
  }
  return(num_vec.back()/denom_vec.back());
}

/*
void
deriv_coeffs_rational_polyomial(const vector<double> &num_coeffs,
                                const vector<double> &denom_coeffs,
                                vector<double> &deriv_num_coeffs,
                                vector<double> &deriv_denom_coeffs){
  for(size_t j = 0; j < 2*denom_coeffs.size(); j++){
    double holding_coeff = 0.0;
    for(size_t k = 0; k <= j; k++){
      if(j-k < denom_coeffs.size())
        holding_coeff += denom_coeffs[j-k]*denom_coeffs[k];
    }
    deriv_denom_coeffs.push_back(holding_coeff);
    
  }
  
  for(size_t j = 0; j < num_coeffs.size() + denom_coeffs.size() - 1; j++){
    double holding_coeff = 0.0;
    for(size_t k = 0; k < = j; k++){
      if(k < num_coeffs.size() - 1 && j-k < denom_coeffs.size())
        holding_coeff += num_coeffs[k+1]*(k+1)*denom_coeffs[j-k];
      if(k < num_coeffs.size() && j-k < denom_coeffs.size() - 1)
        holding_coeff -= num_coeffs[k]*denom_coeffs[j-k+1]*(j-k+1);
    }
    deriv_num_coeffs.push_back(holding_coeff);
  }
}
      

double
pade_upperbound_librarysize(const vector<double> &vals_hist,
                            const double upper_bound, //use pade_upperbound_librarysize(.) to compute
                            const double time_step,
                            const double max_time,
                            const double tolerance,
                            size_t &max_terms){ 
  vector<double> coeffs(max_terms, 0.0);
  for(size_t j = 0; j < max_terms; j++)
    coeffs[j] = vals_hist[j+1]*pow(-1, j+2);
  
  const double values_size = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
  
  if(max_terms % 2 == 1)
    max_terms--; //need max_terms = L+M+1 to be even so that L+M is odd so that we have convergence from below
  //easier to compute derivative using pade
  size_t num_size = max_terms/2;
  size_t denom_size = max_terms - num_size;
  assert(num_size == denom_size);
  vector<double> num_coeffs;
  vector<double> denom_coeffs;
  compute_pade_coeffs(coeffs, num_size, denom_size, num_coeffs, denom_coeffs);
  vector<double> full_denom_coeffs(denom_coeffs);
  full_denom_coeffs.insert(full_denom_coeffs.begin(), 1.0);
  
  vector<double> deriv_num_coeffs;
  vector<double> deriv_denom_coeffs;
  deriv_coeffs_rational_polyomial(num_coeffs, full_denom_coeffs, deriv_num_coeffs, deriv_denom_coeffs);
  vector<double> possible_maxima;
  double t = time_step;
  double current_val  = evaluate_polynomial(deriv_num_coeffs, t);
  double prev_val = evaluate_polynomial(deriv_num_coeffs, 0.0);
  while(t < max_time){
    current_val = evaluate_polynomial(deriv_num_coeffs, t);
    if(current_val*prev_val < 0.0) //detect change
      possible_maxima.push_back(locate_polynomial_zero(deriv_num_coeffs, t, t-time_step, tolerance));
    prev_val = current_val;
    t += time_step;
  }
  possible_maxima.push_back(max_time); //include boundary
  double current_global_max = 0.0;
  double current_global_max_loc = 0.0;
  for(size_t j = 0; j < possible_maxima.size(); j++){
    double test_val = compute_pade_approx_numerator(possible_maxima[j], num_coeffs)/
    compute_pade_approx_denominator(possible_maxima[j], denom_coeffs);
    if(test_val > current_global_max && test_val < upper_bound){
      current_global_max = test_val;
      current_global_max_loc = possible_maxima[j];
    }
  }
  return(current_global_max);
}
*/
                 

double
pade_lowerbound_librarysize_complex(const vector<double> &vals_hist,
                                    const double upper_bound, //use pade_lowerbound_librarysize(.) to compute
                                    const double time_step,
                                    const double max_time,
                                    const double dx,
                                    const double tolerance,
                                    const bool VERBOSE,
                                    size_t max_terms){ 
  double vals_sum = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++)
    vals_sum += i*vals_hist[i];
  
  const double distinct_vals = accumulate(vals_hist.begin(), vals_hist.end(), 0.0); 
  
  vector<double> cf_coeffs;
  vector<double> offset_coeffs;
  cont_frac cf_estimate(cf_coeffs, offset_coeffs, 2, 0);
  
  vector<double> coeffs(max_terms, 0.0);
  for(size_t j = 0; j < max_terms; j++)
    coeffs[j] = vals_hist[j+1]*pow(-1, j+2);
  vector<double> possible_maxima_loc;
  vector<double> possible_maxima;
  if(max_terms % 2 == 0)
    max_terms--; //need max_terms = L+M+1 to be even so that L+M is odd so that we have convergence from below
  while (max_terms > 10){
    if(VERBOSE)
    cf_estimate.compute_cf_coeffs(coeffs, max_terms);
    double t = time_step;
    double prev_deriv = 0.0;
    double current_deriv = cf_estimate.cf_deriv_complex(t, dx, tolerance);
    double current_val = distinct_vals;
    double prev_val = distinct_vals;
    while(t < max_time){
      current_deriv = cf_estimate.cf_deriv_complex(t, dx, tolerance);
      current_val = cf_estimate.cf_approx(t, tolerance)+distinct_vals;
      if(fabs(current_deriv) > vals_sum || current_val > upper_bound){ //if derivative or estimate is not acceptable, choose different order approx
        possible_maxima_loc.clear();  //so that we can know we need to go down in approx
        possible_maxima.clear();
        break; //out of t loop
      }
      else if(current_deriv*prev_deriv < 0.0 && current_val < upper_bound){
        possible_maxima_loc.push_back(cf_estimate.locate_zero_cf_deriv(t, t-time_step,
                                                                       dx, tolerance));
        possible_maxima.push_back(cf_estimate.cf_approx(possible_maxima_loc.back(), tolerance)+distinct_vals);
      }
      else if(current_val < prev_val && prev_val < upper_bound){
        possible_maxima_loc.push_back(t-time_step);
        possible_maxima.push_back(prev_val);
      }
      prev_deriv = current_deriv;
      prev_val = current_val;
      t += time_step;
    }
    if(possible_maxima.size() > 0){
      break; //out of max_terms loop
    }
    else{
      possible_maxima.clear();
      vector<double> no_cf_coeffs;
      vector<double> no_offset_coeffs;
      cf_estimate.set_cf_coeffs(no_cf_coeffs);
      cf_estimate.set_offset_coeffs(no_offset_coeffs);
      max_terms -= 2;
    }
  }
    
  possible_maxima_loc.push_back(max_time); //include boundary
  possible_maxima.push_back(cf_estimate.cf_approx(max_time, tolerance)+distinct_vals);

  if(VERBOSE){
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
    double test_val = possible_maxima[j];
    if(test_val > current_global_max && test_val < upper_bound ){ 
      current_global_max = test_val;
      current_global_max_loc = possible_maxima_loc[j];
    }
  }
  if(VERBOSE)
    cerr << "chosen global max = " << current_global_max << ", loc = " << current_global_max_loc << "\n";
  return(current_global_max);
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
    bool coverage = false;
    bool library_size = false;
    bool smooth_hist = false; 
    double delta = 1e-8;

    
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
    opt_parse.add_opt("max_terms",'x',"maximum number of terms to use for approx",
                      false, max_terms);
    opt_parse.add_opt("delta",'d',"step size for derivative",
                      false, delta);
    opt_parse.add_opt("smooth_hist",'h',"include to smooth histogram",
                      false, smooth_hist);
    opt_parse.add_opt("coverage",'\0', "estimate coverage, default estimates distinct",
                      false, coverage);
    opt_parse.add_opt("library_size",'\0',
                      "estimate upper and lower bound for library size, set max_time to largest time for search for global max",
                      false, library_size);
    

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
    
    const double distinct_vals = static_cast<double>(values.size());
    const size_t vals_sum = accumulate(values.begin(), values.end(), 0);

    const size_t max_value = *std::max_element(values.begin(),values.end());
    vector<double> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i){
      ++vals_hist[static_cast<size_t>(values[i])];
    }
    
    size_t max_approx = 1;
    while(vals_hist[max_approx] > 0)
      max_approx++;
    
    if(smooth_hist){
      size_t smoothing_bandwidth = 4;
      double smoothing_decay_factor = 15.0;
      smooth_histogram(smoothing_bandwidth, smoothing_decay_factor, vals_hist);
    }
    
    vector<double> estimate_vec;
    
    //max_terms selection needs to be dependent on max time
    // heuristics suggest the method is accurate even for large numbers,
    // set (max_time)^max_term < 10^250
    max_terms = std::min(max_approx, max_terms);
    if(VERBOSE){
      cerr << "max term = " << max_terms << "\n";
      cerr << "distinct_vals = " << distinct_vals << "\n";
    }
    
    
    if(coverage){
    //  compute_coverage(vals_hist, max_time, time_step, tolerance,
    //                   defect_error, VERBOSE, estimate_vec, max_terms);
    }
    else if(library_size){
      estimate_vec.push_back(chao_lee_lowerbound_librarysize(vals_hist));
      estimate_vec.push_back(chao87_lowerbound_librarysize(vals_hist));
      double upper_bound = pade_upperbound_librarysize(vals_hist, VERBOSE, max_terms)+distinct_vals;
      estimate_vec.push_back(pade_lowerbound_librarysize_complex(vals_hist, upper_bound, time_step,
                                                                 max_time, delta, tolerance, VERBOSE,
                                                                 max_terms));
      estimate_vec.push_back(upper_bound);

      ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
      *out << "Chao-Lee(1992) lower bound\t" << estimate_vec[0] << endl;
      *out << "Chao(1987) lower bound \t" << estimate_vec[1] << endl;
      *out << "continued fraction lower bound\t" << estimate_vec[2] << endl;
      *out << "continued fraction upper bound\t" << estimate_vec[3] << endl;
      
    }
    else{
      compute_distinct(vals_hist, vals_sum, max_time,
                      time_step, delta, tolerance, VERBOSE,
                       estimate_vec, max_terms);
      double t = 0.0;
      if(estimate_vec.size() > 0){
        ostream* out = (outfile.empty()) ? 
        &std::cout : new std::ofstream(outfile.c_str());
        for(size_t i = 0; i < estimate_vec.size(); i++){
          *out << (t+1.0)*vals_sum << "\t" << estimate_vec[i] << endl;
          t += time_step;
        }
      }
      else{
        cerr << "no acceptable approx \n";
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
