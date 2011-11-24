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


#include <OptionParser.hpp>
#include <GenomicRegion.hpp>
#include <smithlab_utils.hpp>

#include "pade_approximant.hpp"
#include "continued_fraction.hpp"

#include <iostream>
#include <numeric>
#include <vector>

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

using std::greater;

static void
get_counts(const vector<SimpleGenomicRegion> &reads, vector<size_t> &values) {
  size_t count = 1;
  for (size_t i = 1; i < reads.size(); i++)
    if (reads[i] == reads[i-1])
      ++count;
    else {
      values.push_back(count);
      count = 1;
    }
  values.push_back(count);
}

int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    
    double tolerance = 1e-20;
    size_t max_terms = 100;
    double max_time = 1.0;
    double time_step = 1.0;
    double defect_tolerance = 1e-4;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("distinct", "", "<bed-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("time",'s',"maximum time",
		      false, max_time);
    opt_parse.add_opt("time_step",'p',"time step",
                      false, time_step);
    opt_parse.add_opt("terms",'t',"maximum number of terms",
                      false, max_terms);
    opt_parse.add_opt("tol", '\0', "general numerical tolerance",
		      false, tolerance);
    opt_parse.add_opt("def",'\0',"defect tolerance = abs(zero - pole)", 
		      false, defect_tolerance);
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
    
    vector<size_t> values;
    get_counts(read_locations, values);
    
    const size_t vals_sum = accumulate(values.begin(), values.end(), 0);
    
    const size_t max_value = *std::max_element(values.begin(), values.end());
    vector<size_t> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++vals_hist[static_cast<size_t>(values[i])];
    
    max_terms = std::min(max_terms, max_value);
    
    // need max_terms = L+M+1 to be even so that L+M is odd and we get
    // convergence from above
    if(max_terms % 2 == 0)
      max_terms--; 
    vector<double> summands;
    vector<double> coeffs(max_terms, 0.0);
    for(size_t j = 0; j < max_terms; j++)
      coeffs[j] = vals_hist[j+1]*pow(-1, j+2);
    
    size_t denom_size = max_terms/2;
    
    vector<double> numerator_approx;
    vector<double> denominator_approx;
    
    bool DEFECT_FLAG = false;
    size_t max_number_approxs = static_cast<size_t>(round((max_terms-8)/2));
    for(size_t i = 0; i < max_number_approxs; i++) { 
      //detect defect, move up table to more conservative approx
      compute_pade_curve(coeffs, max_time, time_step, defect_tolerance, tolerance,
                         denom_size, VERBOSE, numerator_approx, denominator_approx, 
			 DEFECT_FLAG);
      if (DEFECT_FLAG == false)
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
    
    // if moving up table doesn't work, move down table
    if (DEFECT_FLAG == true) {
      while (coeffs.size() < vals_hist.size()-1) {
        compute_pade_curve(coeffs, max_time, time_step, defect_tolerance, tolerance,
                           denom_size, VERBOSE, numerator_approx, denominator_approx, 
			   DEFECT_FLAG);
        if (DEFECT_FLAG == false)  //exit if no defect
          break;
        coeffs.push_back(vals_hist[coeffs.size()+1]);
        coeffs.push_back(vals_hist[coeffs.size()+1]);
        denom_size++;
        numerator_approx.clear();
        denominator_approx.clear();
        DEFECT_FLAG = false; //reset flag for next iter
      }
    }
    
    // exit if no acceptable approximation found
    if (DEFECT_FLAG) {
      cerr << "no acceptable approximation" << endl; 
      return EXIT_SUCCESS;
    }
    

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    double t = 0.0;
    for(size_t i = 0; i < numerator_approx.size(); i++){
      out << t << "\t" 
	  << (t+1.0)*vals_sum << "\t" 
	  << numerator_approx[i] << "\t"
	  << denominator_approx[i] << "\t" 
	  << t*numerator_approx[i]/denominator_approx[i] + values.size()
	  << endl;
      t += time_step;
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
