/*    poiss_estimation_extrapolation:
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



#include "Poiss_mixture.hpp"

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


#include <fstream>
#include <iomanip>
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

using std::max;
using std::setw;
using std::fabs;
using std::ceil;
using std::greater;
using std::numeric_limits;


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

static inline double
get_max(const double x1, const double x2){
  double larger = (x1 > x2) ? x1 : x2;
  return(larger);
}

static inline double
calculate_neg_AIC(const size_t number_mixtures,
                  const double log_like){
  return(2*log_like - 2*log(2*number_mixtures - 1));
}

static inline double
calculate_neg_BIC(const size_t number_mixtures,
                  const size_t sample_size,
                  const double log_like){
  return(2*log_like - (2*number_mixtures-1)*log(sample_size));
}

static inline size_t
select_number_mixtures(const vector<double> log_like){
  size_t current_max=0;
  double current_max_score = -std::numeric_limits<double>::max();
  for(size_t i = 0; i < log_like.size(); i++){
    double holding_score = calculate_neg_AIC(i+1, log_like[i]);
    if(holding_score >= current_max_score){
      current_max_score = holding_score;
      current_max = i+1;
    }
  }
  return(current_max);
}
   

/*use gsl matrix to do this, i know my code will be slow */
static void
invert_Fisher_info(vector< vector<double> > &Fisher_info){
  gsl_matrix *fish_mat = gsl_matrix_alloc(Fisher_info.size(), 
                                          Fisher_info[0].size());

  for(size_t i = 0; i < Fisher_info.size(); i++){
    for(size_t j = 0; j < Fisher_info[i].size(); j++){
      gsl_matrix_set(fish_mat,i,j, Fisher_info[i][j]);
    }
  }
  gsl_permutation *perm = gsl_permutation_calloc(Fisher_info.size());
  gsl_matrix *inv_fish_mat = gsl_matrix_alloc(Fisher_info.size(),
                                              Fisher_info[0].size());
  int s;
  gsl_linalg_LU_decomp(fish_mat, perm, &s);
  gsl_linalg_LU_invert(fish_mat, perm, inv_fish_mat);
  for(size_t i = 0; i < Fisher_info.size(); i++){
    for(size_t j = 0; j < Fisher_info[i].size(); j++){
      Fisher_info[i][j] = gsl_matrix_get(inv_fish_mat, i, j);
    }
  }
}


int
main(int argc, const char **argv){


  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    size_t max_number_mixtures = 20;
    size_t min_number_mixtures = 1;
    size_t step_size = 1;
    size_t max_iter = 100000;
    double tolerance = 1e-20;
    size_t extrapolation_size = 100000000;
    size_t step_btwn_extra = 1000000;
 

    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("poisson_estimation_extrapolation", "",
			   "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("max_number_mixtures", 'm', 
                      "number of mixtures to take",
                      false, max_number_mixtures);
    opt_parse.add_opt("min_number_mixtures", 'n',
		      "starting number of mixtures",
		      false, min_number_mixtures);
    opt_parse.add_opt("tolerance", 't', "Numerical tolerance",
                     false, tolerance);
    opt_parse.add_opt("step_size",'p',"Steps between mixtures",
		      false, step_size);
    opt_parse.add_opt("extrapolation_size",'e',"Length of extrapolation",
		      false, extrapolation_size);
    opt_parse.add_opt("step_btwn_extra",'b',
		      "Length between extrapolation points",
		      false,step_btwn_extra);


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

    const size_t max_value = *std::max_element(values.begin(),values.end());
    vector<size_t> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i){
      ++vals_hist[static_cast<size_t>(values[i])];
    }
    vector<double> score_vec;
    vector< vector<double> > lambdas_vec;
    vector< vector<double> > mixings_vec;
    vector< vector< vector<double> > > Fisher_info_vec;
    vector< vector<double> > expected_MN_vec;

    for(size_t k = min_number_mixtures; 
	k <= max_number_mixtures; k += step_size){
      size_t number_mixtures = k;
      vector<ZTP> distros;
      vector<double> mixing(number_mixtures, 
                          1/static_cast<double>(number_mixtures));
    
      sort(values.begin(), values.end());
      double step_size = floor(static_cast<double>(values.size())
                                 /static_cast<double>(number_mixtures));
      for(size_t j = 0; j < number_mixtures; j++){
        double lambda_hat = accumulate(values.begin()+j*step_size,
                              values.begin()+(j+1)*step_size, 
                              0.0)/step_size;

        if(j != 0){
          ZTP last_distro = distros.back();
          lambda_hat = get_max(lambda_hat, last_distro.get_lambda()) 
                         + 0.01;
        }
        else
          lambda_hat = lambda_hat+0.01;

 
        ZTP holding_distro(lambda_hat);
        distros.push_back(holding_distro);
      }
      vector< vector<double> > Fish_info(2*k-1, vector<double>(2*k-1, 0.0));
      ZTP_mixture ZTP_mix(distros, mixing, Fish_info);

      double score = ZTP_mix.EM_mix_resolve(vals_hist, tolerance, max_iter);
      distros = ZTP_mix.get_distros();
      vector<double> lambdas;
      for(size_t j = 0; j < distros.size(); j++)
        lambdas.push_back(distros[j].get_lambda());

      lambdas_vec.push_back(lambdas);
      mixings_vec.push_back(ZTP_mix.get_mixing());
      score_vec.push_back(score);
      Fisher_info_vec.push_back(ZTP_mix.get_Fish_info());

    vector<double> expect_MN;
    for(size_t i = step_btwn_extra; i <= extrapolation_size;
	i += step_btwn_extra){
      expect_MN.push_back(ZTP_mix.expected_inverse_sum(values_size, i));
    }
    expected_MN_vec.push_back(expect_MN);
			   
    }

    size_t opt_num = select_number_mixtures(score_vec);
    size_t opt_num_mixs = step_size*(opt_num-1) + min_number_mixtures;

    vector<double> lambdas = lambdas_vec[opt_num-1];
    vector<double> mixing = mixings_vec[opt_num-1];
    vector< vector<double> > Fisher_info=Fisher_info_vec[opt_num-1];
    vector<double> expected_MN = expected_MN_vec[opt_num-1];
    invert_Fisher_info(Fisher_info);

    cerr << opt_num << "\n";

    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    size_t time_step = step_btwn_extra;
    for(size_t i = 0; i < expected_MN.size(); i++){
      *out << time_step << "\t" << expected_MN[i] << endl;
      time_step += step_btwn_extra;
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
