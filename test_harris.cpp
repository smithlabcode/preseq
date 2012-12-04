/*    test_harris:
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

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <smithlab_os.hpp>

#include "library_size_estimates.hpp"
#include "newtons_method.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::fixed;

/*
// want #reads in hist_in = reads in hist_out
static void
renormalize_hist(const vector<double> &hist_in,
                 vector<double> &hist_out) {
  double out_vals_sum = 0.0;
  for(size_t i = 0; i < hist_out.size(); i++)
    out_vals_sum += i*hist_out[i];
  double in_vals_sum = 0.0;
  for(size_t i = 0; i < hist_in.size(); i++)
    in_vals_sum += i*hist_in[i];
  for(size_t  i =0; i < hist_out.size(); i++)
    hist_out[i] = hist_out[i]*in_vals_sum/out_vals_sum;
}
*/

void
randomize_coefficients_compare_bounds(const bool VERBOSE,
				      const vector<double> &counts_hist, 
				      const double tolerance,
				      const size_t max_iter, 
				      const size_t max_terms,
				      const size_t bootstraps,
				      const size_t rand_indx,
				      vector<double> &chao_bounds,
				      vector<double> &harris_bounds){

  chao_bounds.clear();
  harris_bounds.clear();

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand()); 

  while(chao_bounds.size() < bootstraps){

    vector<double> rand_counts_hist(counts_hist);

    rand_counts_hist[rand_indx] = gsl_ran_poisson(rng, counts_hist[rand_indx]);

    //    renormalize_hist(counts_hist, rand_counts_hist);

    const double harris_lower_bound = my_harris(VERBOSE, rand_counts_hist, tolerance,
						max_iter, max_terms);

    const double chao_lower_bound = chao87_lowerbound_librarysize(rand_counts_hist);

    if(harris_lower_bound > 0.0 && chao_lower_bound > 0.0){
      harris_bounds.push_back(harris_lower_bound);
      chao_bounds.push_back(chao_lower_bound);
      cerr << "converge, iter " << harris_bounds.size() << endl;
    }
  }
}

void
compute_mediandiff_CI(const vector<double> &chao_bounds,
		    const vector<double> &harris_bounds,
		    const double alpha, double &median_diff,
		    double &upper_CI, double &lower_CI){
  vector<double> bounds_diff(harris_bounds);
  assert(bounds_diff.size() == chao_bounds.size());
  for(size_t i = 0; i < bounds_diff.size(); i++)
    bounds_diff[i] -= chao_bounds[i];

  sort(bounds_diff.begin(), bounds_diff.end());
  median_diff = gsl_stats_median_from_sorted_data(&bounds_diff[0], 1, bounds_diff.size());
  lower_CI = gsl_stats_quantile_from_sorted_data(&bounds_diff[0], 1, 
						 bounds_diff.size(), alpha);
  upper_CI = gsl_stats_quantile_from_sorted_data(&bounds_diff[0], 1, 
						 bounds_diff.size(), 1.0 - alpha);
}

int
main(const int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;
    
    size_t orig_max_terms = 5;
    double lib_size = 3.0e9;
    double sampled_reads = 1e7;
    double alpha = 1.0;
    double tolerance = 1e-8;
    size_t max_iter = 1000;
    size_t hist_max_terms = 1000;
    size_t rand_indx = 1;
    size_t bootstraps = 1000;
    double CI = 0.95;

    
    /* FLAGS */
    bool VERBOSE = false;
    //	bool SMOOTH_HISTOGRAM = false;	
    
    /**************** GET COMMAND LINE ARGUMENTS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "",
			   "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "yield output file (default: stdout)",
		      false , outfile);
    opt_parse.add_opt("max_terms",'x', "max terms",
		      false, orig_max_terms);
    opt_parse.add_opt("hist_max_terms",'h',"max terms in histogram",
		      false, hist_max_terms);
    opt_parse.add_opt("lib_size",'l', "library size",
		      false, lib_size);
    opt_parse.add_opt("sampled_reads",'s',"number of sampled reads",
		      false, sampled_reads);
    opt_parse.add_opt("alpha",'a',"alpha for NegBin dist",
		      false, alpha);
    opt_parse.add_opt("tol",'t',"numerical tolerance",
		      false, tolerance);
    opt_parse.add_opt("max_iter",'i',"maximum # iterations",
		      false, max_iter);
    opt_parse.add_opt("rand_indx",'r',"index to randomize",
		      false, rand_indx);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps to perform",
		      false, bootstraps);
    opt_parse.add_opt("CI",'c', "Confidence level",
		      false, CI);
    //	opt_parse.add_opt("terms",'t',"maximum number of terms", false, 
    //	     orig_max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false, VERBOSE);
    
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
    /******************************************************************/
    
    const double alpha_inv = 1/alpha;
    // BUILD THE HISTOGRAM
    double mu = sampled_reads/lib_size; 
    vector<double> counts_hist;
    counts_hist.push_back(0.0);
    for(size_t i = 1; i < hist_max_terms; i++){
      /*  counts_hist.push_back(exp(log(lib_size)  - mu 
				+i*log(mu) - gsl_sf_lnfact(i)));
      */
      counts_hist.push_back(exp(log(lib_size)
			       + gsl_sf_lngamma(i + alpha_inv) - gsl_sf_lnfact(i)
			       - gsl_sf_lngamma(alpha_inv) + i*log(alpha)
			       + i*log(mu) 
				- (i + alpha_inv)*log(1.0 + alpha*mu)));	  
	  
      if(counts_hist.back() < tolerance)
	break;
    }
    vector<double> measure_moments;
  // mu_r = (r + 1)! n_{r+1} / n_1
    for(size_t i = 0; i < std::max(orig_max_terms, counts_hist.size() - 1); i++)
      measure_moments.push_back(exp(gsl_sf_lngamma(i + 3)
				    + log(counts_hist[i + 2])
				    - log(counts_hist[1])));
    // if (VERBOSE) {
      cerr << "LIBRARY_SIZE = " << lib_size << endl;
      cerr << "MU = " << mu << "\t" << "ALPHA = " << alpha << endl;
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << "MOMENTS" << endl;
      for(size_t i = 0; i < measure_moments.size(); i++)
	cerr << measure_moments[i] << endl;
      cerr << "OBSERVED_DISTINCT = " << accumulate(counts_hist.begin(), counts_hist.end(), 0.0) << endl;
      //   }
    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////

    vector<double> chao_bounds, harris_bounds;
    randomize_coefficients_compare_bounds(VERBOSE, counts_hist, tolerance,
					  max_iter, orig_max_terms,
					  bootstraps, rand_indx,
					  chao_bounds, harris_bounds);

    double median_diff, upper_CI, lower_CI;
    compute_mediandiff_CI(chao_bounds, harris_bounds, 1.0 - CI, median_diff,
			upper_CI, lower_CI);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << "median_diff" << '\t' << "lower" << CI << "CI" 
	<< '\t' << "upper" << CI << "CI" << endl;

    out << median_diff << '\t' << lower_CI << '\t'
	<< upper_CI << endl;




    
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////

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
