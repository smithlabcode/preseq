/*    negbin_estimation_extrapolation:
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



#include "NBD_mixture.hpp"

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

static inline double
d_MN_d_real_mixing(const double current_mu,
		   const double current_alpha,
		   const size_t N){
  return(exp(log(N)-log(current_mu)+
	     log(1-exp(-log(1+current_mu*current_alpha)/current_alpha))));
}

static inline double
d_MN_d_current_mu(const double current_mu,
		  const double current_alpha,
		  const double real_mixing,
		  const size_t N){
  const double plus_alpha_mu = 1 + current_alpha*current_mu;
  return(exp(log(real_mixing)+log(N)-log(current_mu)
	     +log(current_alpha+1)
	     -((current_alpha+1)/current_alpha)*log(plus_alpha_mu))
	 -exp(log(real_mixing)+log(N)-2*log(current_mu)
	      +log(1-exp(-log(plus_alpha_mu)/current_alpha))));
}

static inline double
d_MN_d_alpha(const double current_mu,
	     const double alpha,
	     const double real_mixing,
	     const size_t N){
  const double plus_alpha_mu = 1 + alpha*current_mu;
  return(-exp(log(real_mixing)+log(N)-log(current_mu)
	      -log(plus_alpha_mu)/alpha-2*log(alpha))*log(plus_alpha_mu)
	 +exp(log(real_mixing)+log(N)-log(current_mu)
	      -log(plus_alpha_mu)/alpha-log(alpha)
	      -log(plus_alpha_mu)));
}

static inline double
d_current_mu_d_mu_same_indx(const double mu,
			    const double alpha,
			    const double theta,
			    const double mean,
			    const size_t N,
			    const size_t sample_size){
  const double denom = 1 - exp(-log(1+alpha*mu)/alpha);
  return(exp(log(N)-log(sample_size) - log(mean))
	 -exp(log(N)-log(sample_size)-2*log(mean))*
	 (theta/denom - exp(log(theta)+log(mu)-2*log(denom)
			    -(alpha+1)*log(1+alpha*mu)/alpha)));
}

static inline double
d_current_mu_dmu_diff_indx(const double mu_num,
			   const double mu_denom,
			   const double alpha_denom, 
			   const double theta_denom,
			   const double mean, 
			   const size_t N,
			   const size_t sample_size){
  const double denom = 1 - exp(-log(1+alpha_denom*mu_denom)/alpha_denom);
  return(-exp(log(N)-log(sample_size)+log(mu_num)-2*log(mean))*
	 (theta_denom/denom - 
	  exp(log(theta_denom)+log(mu_denom)-2*log(denom)
	      -(alpha_denom+1)*log(1+alpha_denom*mu_denom)/alpha_denom)));
}

static inline double
d_current_mu_dtheta_all_indx(const double mu_mu,
			     const double mu_theta,
			     const double mu_last,
			     const double alpha_theta,
			     const double alpha_last,
			     const double mean,
			     const size_t N,
			     const size_t sample_size){
  const double denom_theta = 1-exp(-log(1+alpha_theta*mu_theta)/alpha_theta);
  const double denom_last = 1-exp(-log(1+alpha_last*mu_last)/alpha_last);
  return(exp(log(N)-log(sample_size)+log(mu_mu)-2*log(mean))*
	 (-mu_theta/denom_theta + mu_last/denom_last));
}


static void
resamplevals(const gsl_rng *rng,
	     const vector<size_t> &orig,
	     const vector<size_t> &vals,
	     vector<size_t> &sample){
  vector<size_t> indx_sample(sample.size(), 0);
  gsl_ran_sample(rng, (size_t *)&indx_sample.front(), sample.size(), 
		 (size_t *)&orig.front(), orig.size(), sizeof(size_t));
  
  for(size_t i = 0; i < sample.size(); i++)
    sample[i] = vals[indx_sample[i]];
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
    size_t bootstraps = 0;
    bool sampling_error = false;
 

    
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
    opt_parse.add_opt("step_btwn_extra",'w',
		      "Length between extrapolation points",
		      false,step_btwn_extra);
    opt_parse.add_opt("bootstraps",'b',
		      "Number of bootstraps for CI, exclude for normal CI",
		      false, bootstraps);
    opt_parse.add_opt("sampling_error",'s',
		      "Include sampling variance in CI",
		      false, sampling_error);


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
    std::ifstream in(input_file_name.c_str());
    if (!in) 
      throw BEDFileException("cannot open input file " + input_file_name);
    static const size_t buffer_size = 10000; // Magic!
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    double real_score = atof(buffer);
    while (!in.eof()){
      int holding_val = 0;
      char buffer[buffer_size];
      in.getline(buffer, buffer_size);
      holding_val = atoi(buffer);
      if(!(holding_val >= 0)){
        cerr << "Invalid Input.\n";
      }
      assert(holding_val >= 0);
      values.push_back(holding_val);
    }
    in.close();
    values.pop_back();


    const size_t values_size = values.size();
    const size_t max_value = *std::max_element(values.begin(),values.end());
    vector<size_t> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i){
      ++vals_hist[static_cast<size_t>(values[i])];
    }
    vector<double> score_vec;
    vector< vector<double> > mus_vec;
    vector< vector<double> > mixings_vec;
    vector< vector< vector<double> > > Fisher_info_vec;
    vector< vector<double> > expected_MN_vec;

    for(size_t k = min_number_mixtures; 
	k <= max_number_mixtures; k += step_size){
      size_t number_mixtures = k;
      vector<ZTNBD> distros;
      vector<double> mixing(number_mixtures, 
                          1/static_cast<double>(number_mixtures));
    
      sort(values.begin(), values.end());
      double step_size = floor(static_cast<double>(values.size())
                                 /static_cast<double>(number_mixtures));
      for(size_t j = 0; j < number_mixtures; j++){
        double mu_hat = accumulate(values.begin()+j*step_size,
                              values.begin()+(j+1)*step_size, 
                              0.0)/step_size;

        if(j != 0){
          ZTNBD last_distro = distros.back();
          mu_hat = get_max(mu_hat, last_distro.get_mu()) 
                         + 0.01;
        }
        else
          mu_hat = mu_hat+0.01;

 
        ZTNBD holding_distro(mu_hat, 1);
        distros.push_back(holding_distro);
      }
      vector< vector<double> > Fish_info(2*k-1, vector<double>(2*k-1, 0.0));
      ZTNBD_mixture ZTNBD_mix(distros, mixing, Fish_info);

      double score = 
	ZTNBD_mix.EM_resolve_mix_add_zeros(tolerance, max_iter, vals_hist);

      distros = ZTNBD_mix.get_distros();
      vector<double> mus;
      for(size_t j = 0; j < distros.size(); j++)
        mus.push_back(distros[j].get_mu());

      vector<double> mix = ZTNBD_mix.get_mixing();
      vector< vector<double> > Fish = ZTNBD_mix.get_Fish_info();

      for(size_t j = 0; j < distros.size(); j++){
	cerr << mix[j] << ", " << distros[j].get_mu() 
	     << ", " << distros[j].get_alpha() << "\n";
      }
      cerr << "\n";
      for(size_t j = 0; j < Fish.size(); j++){
	for(size_t l = 0; l < Fish[j].size(); l++){
	  cerr << Fish[j][l] << "\t";
	}
	cerr << "\n";
      }

      mus_vec.push_back(mus);
      mixings_vec.push_back(ZTNBD_mix.get_mixing());
      score_vec.push_back(score);
      Fisher_info_vec.push_back(ZTNBD_mix.get_Fish_info());


    vector<double> expect_MN;
    for(size_t i = step_btwn_extra; i <= extrapolation_size;
	i += step_btwn_extra){
      expect_MN.push_back(ZTNBD_mix.expected_inverse_sum(values_size, i));
    }
    expected_MN_vec.push_back(expect_MN);
			   
    }

    size_t opt_num = select_number_mixtures(score_vec);
    size_t opt_num_mixs = step_size*(opt_num-1) + min_number_mixtures;

    vector<double> mus = mus_vec[opt_num-1];
    vector<double> mixing = mixings_vec[opt_num-1];
    vector< vector<double> > Fisher_info=Fisher_info_vec[opt_num-1];
    vector<double> expected_MN = expected_MN_vec[opt_num-1];


    /*
    if(bootstraps == 0){
      invert_Fisher_info(Fisher_info);


      vector<double> current_lambdas(lambdas);
      vector<double> real_mixing(mixing);
      calculate_real_mixing(lambdas, mixing, real_mixing);

      if(sampling_error == false){
	ostream* out = (outfile.empty()) ? 
	  &std::cout : new std::ofstream(outfile.c_str());
	size_t time_step = step_btwn_extra;
	for(size_t i = 0; i < expected_MN.size(); i++){
	  calculate_current_lambdas(lambdas, real_mixing, mixing, values_size,
				    time_step, current_lambdas);


	  double log_var_fit = 
	    compute_log_var_from_fit(Fisher_info, lambdas,mixing,
				     current_lambdas,real_mixing,
				     time_step, values_size);
	  *out << time_step << "\t" << expected_MN[i] << "\t"
	       << expected_MN[i]-1.959964*exp(0.5*log_var_fit) << "\t"
	       << expected_MN[i]+1.959964*exp(0.5*log_var_fit) << endl;
	  time_step += step_btwn_extra;
	}
      }
   
      else{
	ostream* out = (outfile.empty()) ? 
	  &std::cout : new std::ofstream(outfile.c_str());
	size_t time_step = step_btwn_extra;
	vector<double> current_mixing(mixing);
	double var_samp = 0.0;

	for(size_t i = 0; i < expected_MN.size(); i++){
	  calculate_current_lambdas(lambdas, real_mixing, mixing, values_size,
				    time_step, current_lambdas);

	  calculate_current_mixings(real_mixing, current_lambdas, 
				    time_step, current_mixing);
	  if(i % 10 == 0){
	  var_samp = 
	    exp(compute_log_var_from_sampling(current_lambdas,current_mixing,
					      real_mixing,
					      time_step+10*step_btwn_extra,
					      tolerance));
	  }
	  double var_fit = 
	    exp(compute_log_var_from_fit(Fisher_info, lambdas,mixing,
				     current_lambdas,real_mixing,
					 time_step, values_size));
	  *out << time_step << "\t" << expected_MN[i] << "\t"
	       << expected_MN[i]-1.959964*exp(0.5*(log(var_fit+var_samp))) 
	       << "\t"
	       << expected_MN[i]+1.959964*exp(0.5*(log(var_fit+var_samp))) 
	       << endl;
	  time_step += step_btwn_extra;
	}
      }
    }
    else{

      vector<size_t> orig(values.size());
      for (size_t i = 0; i < orig.size(); ++i)
	orig[i] = i;

      gsl_rng_env_setup();
      gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
      srand(time(0) + getpid());
      gsl_rng_set(rng, rand());
 
      vector< vector<double> > bootstrap_MN;
      size_t number_mixtures = opt_num_mixs;
      for(size_t i = 0; i < bootstraps; i++){

	vector<size_t> sample(values.size(), 0);
	resamplevals(rng,orig,values,sample);
	const size_t max_val = *std::max_element(sample.begin(),sample.end());
	vector<size_t> sample_hist(max_val + 1, 0.0);
	for (size_t j = 0; j < sample.size(); ++j){
	  ++sample_hist[static_cast<size_t>(sample[j])];
	}
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
	vector< vector<double> > Fish_info(2*number_mixtures-1, 
					   vector<double>(2*number_mixtures-1, 
							  0.0));
	ZTP_mixture ZTP_mix(distros, mixing, Fish_info);

	double score = ZTP_mix.EM_mix_resolve(vals_hist, tolerance, max_iter);
	vector<double> expect_MN;
	for(size_t i = step_btwn_extra; i <= extrapolation_size;
	    i += step_btwn_extra){
	  expect_MN.push_back(ZTP_mix.expected_inverse_sum(values_size, i));
	}
	bootstrap_MN.push_back(expect_MN);			  
      }
      //sort bootstrapped data
      for(size_t i = 0; i < bootstrap_MN[0].size(); i++){
	vector<double> holding_vec;
	for(size_t  j= 0; j < bootstrap_MN.size(); j++){
	  holding_vec.push_back(bootstrap_MN[j][i]);
	}
	sort(holding_vec.begin(), holding_vec.end());
	for(size_t j = 0; j < bootstrap_MN.size(); j++){
	  bootstrap_MN[j][i] = holding_vec[j];
	}
      }

      for(size_t i = 0; i < bootstrap_MN[0].size(); i++){
	vector<double> holding_vec;
	for(size_t  j= 0; j < bootstrap_MN.size(); j++){
	  holding_vec.push_back(bootstrap_MN[j][i]);
	}
	sort(holding_vec.begin(), holding_vec.end());
	for(size_t j = 0; j < bootstrap_MN.size(); j++){
	  bootstrap_MN[j][i] = holding_vec[j];
	}
      }
      ostream* out = (outfile.empty()) ? 
	&std::cout : new std::ofstream(outfile.c_str());
      for(size_t i = 0; i < bootstrap_MN.size(); i++){
	for(size_t j = 0; j < bootstrap_MN[i].size(); j++){
	  *out << bootstrap_MN[i][j] << "\t";
	}
	*out << endl;
      }  
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
