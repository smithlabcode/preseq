/*    distinct:
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




int
main(int argc, const char **argv){


  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    double tolerance = 1e-20;
    size_t max_terms = 100;
    double max_time = 1.0;
    double time_step = 1.0;
    double defect_error = 1e-4;
 

    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("distinct", "",
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
    opt_parse.add_opt("max_terms",'x',"maximum number of terms",
                      false, max_terms);
    opt_parse.add_opt("defect_error",'e',"allowable error between zero and pole in approximation",
                      false, defect_error);
    




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
    
    size_t values_size = values.size();
    size_t vals_sum = accumulate(values.begin(), values.end(), 0);

    const size_t max_value = *std::max_element(values.begin(),values.end());
    vector<size_t> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i){
      ++vals_hist[static_cast<size_t>(values[i])];
    }
    
    max_terms = std::min(max_terms, max_value);
    
    if(max_terms % 2 == 0)
      max_terms--; //need max_terms = L+M+1 to be even so that L+M is odd and we get convergence from above
    vector<double> summands;
    vector<double> coeffs(max_terms, 0.0);
    for(size_t j = 0; j < max_terms; j++)
      coeffs[j] = vals_hist[j+1]*pow(-1, j+2);
    
    size_t denom_size = max_terms/2; //denom_size = M
    
    vector<double> numerator_approx;
    vector<double> denominator_approx;
    
    bool defect_flag = false;
    size_t max_number_approxs = static_cast<size_t>(round((max_terms-8)/2));
    for(size_t i = 0; i < max_number_approxs; i++){ //detect defect, move up table to more conservative approx
      compute_pade_curve(coeffs, max_time, time_step, defect_error, tolerance,
                         denom_size, numerator_approx, denominator_approx, defect_flag);
      if(defect_flag == false)
        break;  //no allowable defect detected, accept approx
      coeffs.pop_back();
      coeffs.pop_back();
      denom_size--;
      numerator_approx.clear();
      denominator_approx.clear();
      cerr << i+1 << "\n";
    }
    cerr << "\n";
    
    coeffs.clear();
    for(size_t j = 0; j < max_terms+2; j++)
      coeffs.push_back(vals_hist[j+1]);
    denom_size = coeffs.size()/2 - 1; 
    
    //if moving up table doesn't work, move down table
    if(defect_flag == true){
      while(coeffs.size() < vals_hist.size()-1){
        compute_pade_curve(coeffs, max_time, time_step, defect_error, tolerance,
                           denom_size, numerator_approx, denominator_approx, defect_flag);
        if(defect_flag == false)  //exit if no defect
          break;
        coeffs.push_back(vals_hist[coeffs.size()+1]);
        coeffs.push_back(vals_hist[coeffs.size()+1]);
        denom_size++;
        numerator_approx.clear();
        denominator_approx.clear();
        defect_flag = false; //reset flag for next iter
      }
    }
    
    if(defect_error == true){
      cerr << "no allowable approximation" << endl; //exit if no acceptable approximation found
      return EXIT_SUCCESS;
    }
    
    double t = 0.0;
    ostream* out = (outfile.empty()) ? 
    &std::cout : new std::ofstream(outfile.c_str());
    for(size_t i = 0; i < numerator_approx.size(); i++){
      *out << t << "\t" << (t+1.0)*vals_sum << "\t" << numerator_approx[i] << "\t"
      << denominator_approx[i] << "\t" << t*numerator_approx[i]/denominator_approx[i] + values_size
      << endl;
      t += time_step;
    }
      
    
    
    /*
    size_t numer_size = max_terms - denom_size;  //numer_size = L+1
    
    
    vector<double> denom_vec;
    vector<double> num_vec;
    
    cerr << "numer size = " << numer_size << ", denom size = " << denom_size << "\n";
    
    compute_pade_coeffs(coeffs, numer_size, denom_size, num_vec, denom_vec); 
    
    vector<double> numerator_approx;
    vector<double> denominator_approx;
    
    double t = 0.0;
    
    double prev_denom_val = 1.0;
    double current_denom_val = 1.0;
    double zero_location = 0.0; 
    
    while(t <= max_time){
      numerator_approx.push_back(compute_pade_approx_numerator(t, num_vec));
      denominator_approx.push_back(compute_pade_approx_denominator(t, denom_vec));
    
    
    
    ostream* out = (outfile.empty()) ? 
    &std::cout : new std::ofstream(outfile.c_str());



    while(t <= max_time){
      *out << t << "\t" << (t+1.0)*vals_sum << "\t" << compute_pade_approx_numerator(t, num_vec)
      << "\t" << compute_pade_approx_denominator(t, denom_vec) << "\t"
      << t*compute_pade_approx_numerator(t, num_vec)/compute_pade_approx_denominator(t, denom_vec) + values_size
      << endl;
      current_denom_val = compute_pade_approx_denominator(t, denom_vec);
      if(current_denom_val*prev_denom_val < 0){
        vector<double> denom_coeffs(denom_vec); 
        denom_coeffs.insert(denom_coeffs.begin(), 1.0);
        zero_location = locate_polynomial_zero(denom_coeffs, t-time_step, t, tolerance);
        cerr << "zero located, lower limit = " << t-time_step << ", value = " << compute_pade_approx_denominator(t-time_step, denom_vec)
        << ", upper limit = " << t << ", value = " << compute_pade_approx_denominator(t, denom_vec) 
        << ", location = " << zero_location << "\n";
      }
      prev_denom_val = current_denom_val;
      t += time_step;
    }
      */

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
