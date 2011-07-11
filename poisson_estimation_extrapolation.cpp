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
      current_max = i;
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

void
calculate_real_mixing(const vector<double> &lambdas,
		      const vector<double> &mixing,
		      vector<double> &real_mixing){
  double mean  = 0;
  for(size_t i = 0; i < lambdas.size(); i++)
    mean += mixing[i]*lambdas[i]/(1-exp(-lambdas[i]));
  for(size_t i = 0; i < real_mixing.size(); i++)
    real_mixing[i] = mixing[i]*lambdas[i]/((1-exp(-lambdas[i]))*mean);
}


static inline double
d_MN_d_current_lambda(const double current_lambda,
		      const double real_mixing,
		      const size_t N){
  return(exp(log(N) + log(real_mixing) - current_lambda 
	     - log(current_lambda))
	 -exp(log(N) + log(real_mixing) + log(1-exp(-current_lambda))
	      -2*log(current_lambda)));
}

static inline double
d_MN_d_real_mixing(const double current_lambda,
		   const size_t N){
  return(exp(log(N)+log(1-exp(-current_lambda)) - log(current_lambda)));
}

static inline double
d_current_lambda_d_lambda_same_indx(const double lambda,
				    const double mixing,
				    const double mean,
				    const size_t N,
				    const size_t values_size){
  return(exp(log(N)-log(values_size)-log(mean))
	 -exp(log(N)-log(values_size)-2*log(mean)+log(mixing)
	      -log(1-exp(-lambda)))
	 +exp(log(N)-log(values_size)-2*log(mean)+log(mixing)
	      +log(lambda)-lambda-2*log(1-exp(-lambda))));
}

static inline double
d_current_lambda_d_lambda_dif_indx(const double lambda_deriv,
				   const double lambda_num,
				   const double mixing,
				   const double mean,
				   const size_t N,
				   const size_t values_size){
  return(exp(log(N)-log(values_size) + log(lambda_num) -2*log(mean)
	     +log(mixing)+log(lambda_deriv)-lambda_deriv
	     -2*log(1-exp(-lambda_deriv)))
	 -exp(log(N)-log(values_size)+log(lambda_num)-2*log(mean)
	      +log(mixing)-log(1-exp(-lambda_deriv))));
}

static inline double
d_current_lambda_d_mixing(const double lambda_d_indx,
			  const double lambda_last,
			  const double lambda_current,
			  const double mean,
			  const size_t N,
			  const size_t values_size){
  return(-exp(log(N)-log(values_size)+log(lambda_current)-2*log(mean)
	     +log(lambda_d_indx)-log(1-exp(-lambda_d_indx)))
	 +exp(log(N)-log(values_size)+log(lambda_current)-2*log(mean)
	      +log(lambda_last)-log(1-exp(-lambda_last))));
}

static inline double
d_real_mix_d_mix_same_indx(const double lambda,
			   const double lambda_last,
			   const double mixing,
			   const double mean){
  return(exp(log(lambda)-log(mean)-log(1-exp(-lambda)))
	 -exp(log(mixing)+2*log(lambda)-2*log(1-exp(-lambda))-2*log(mean))
	 +exp(log(mixing)+log(lambda)-log(1-exp(-lambda))
	      -2*log(mean)+log(lambda_last)-log(1-exp(-lambda_last))));
}

static inline double
d_real_mix_d_mix_dif_indx(const double lambda_real,
			  const double lambda_deriv,
			  const double lambda_last,
			  const double mixing_real,
			  const double mean){
  return(-exp(log(mixing_real)+log(lambda_real)-log(1-exp(-lambda_real))
	      -2*log(mean)+log(lambda_deriv)-log(1-exp(-lambda_deriv)))
	 +exp(log(mixing_real)+log(lambda_real)-log(1-exp(-lambda_real))
	      -2*log(mean)+log(lambda_last)-log(1-exp(-lambda_last))));
}

static inline double
d_real_mix_d_lambda_same_indx(const double lambda,
			      const double mixing,
			      const double mean){
  return(exp(log(mixing)-log(1-exp(-lambda))-log(mean))
	 -exp(log(mixing)+log(lambda)-lambda-log(mean)-2*log(1-exp(-lambda)))
	 -exp(2*log(mixing)+log(lambda)-2*log(1-exp(-lambda))-2*log(mean))
	 +exp(2*log(mixing)+2*log(lambda)-lambda-2*log(mean)
	      -3*log(1-exp(-lambda))));
}

static inline double
d_real_mix_d_lambda_dif_indx(const double lambda_real,
			     const double mixing_real,
			     const double lambda_deriv,
			     const double mixing_deriv,
			     const double mean){
  return(exp(log(mixing_real)+log(lambda_real)-2*log(mean)
	     -log(1-exp(-lambda_real))+log(mixing_deriv)
	     +log(lambda_deriv)-lambda_deriv
	     -2*log(1-exp(-lambda_deriv)))
	 -exp(log(mixing_real)+log(lambda_real)-2*log(mean)
	      -log(1-exp(-lambda_real))+log(mixing_deriv)
	      -log(1-exp(-lambda_deriv))));
}

static inline double
d_last_real_mix_d_mix(const double lambda,
		      const double last_lambda,
		      const double last_mixing,
		      const double mean){
  return(-exp(log(last_lambda)-log(mean)-log(1-exp(-last_lambda)))
	 -exp(log(last_mixing)+log(last_lambda)-log(1-exp(-last_lambda))
	      -2*log(mean)+log(lambda)-log(1-exp(-lambda)))
	 +exp(log(last_mixing)+2*log(last_lambda)-2*log(mean)
	      -2*log(1-exp(-last_lambda))));
}



static void
compute_deriv_MN_vec(const vector<double> &lambdas,
		     const vector<double> &current_lambdas,
		     const vector<double> &mixing,
		     const vector<double> &real_mixing,
		     const size_t N,
		     const size_t values_size,
		     vector<double> &MN_deriv){
  vector<double> log_mean_vec;
  for(size_t i = 0; i < lambdas.size(); i++){
    log_mean_vec.push_back(log(mixing[i])+log(lambdas[i])
			   -log(1-exp(-lambdas[i])));
  }
  const double mean = exp(log_sum_log_vec(log_mean_vec,
					  log_mean_vec.size()));

  for(size_t i = 0; i < real_mixing.size()-1; i++){
    MN_deriv[i] = 
      d_MN_d_current_lambda(current_lambdas[i],real_mixing[i],N)
      *d_current_lambda_d_mixing(lambdas[i], lambdas.back(),
				 current_lambdas[i], mean, N, values_size);
    MN_deriv[i] += d_MN_d_real_mixing(current_lambdas[i], N)
      *d_real_mix_d_mix_same_indx(lambdas[i], lambdas.back(), 
				  mixing[i], mean);
    for(size_t j = 0; j < current_lambdas.size(); j++){
      if(j != i){
	MN_deriv[i] += 
	  d_MN_d_current_lambda(current_lambdas[j],real_mixing[j], N)
	  *d_current_lambda_d_mixing(lambdas[j], lambdas.back(),
				     current_lambdas[i], mean, N, values_size);
	if(j != current_lambdas.size()-1){
	  MN_deriv[i] += 
	    d_MN_d_real_mixing(current_lambdas[j], N)
	    *d_real_mix_d_mix_dif_indx(lambdas[j], lambdas[i],
				       lambdas.back(), mixing[j], mean);
	}
	else{
	  MN_deriv[i] += 
	    d_MN_d_real_mixing(current_lambdas[j], N)
	    *d_last_real_mix_d_mix(lambdas[i], lambdas.back(), 
				   mixing.back(), mean);
	}
      }
    }
  }
  for(size_t i = 0; i < current_lambdas.size(); i++){
    size_t indx = i+real_mixing.size()-1;
    MN_deriv[indx] =
      d_MN_d_current_lambda(current_lambdas[i], real_mixing[i], N)
      *d_current_lambda_d_lambda_same_indx(lambdas[i], mixing[i], mean,
					   N, values_size);
    MN_deriv[indx] +=
      d_MN_d_real_mixing(current_lambdas[i], N)
      *d_real_mix_d_lambda_same_indx(lambdas[i], mixing[i], mean);
    for(size_t j = 0; j < current_lambdas.size(); j++){
      if(j != i){
	MN_deriv[indx] += 
	  d_MN_d_current_lambda(current_lambdas[j], real_mixing[j], N)
	  *d_current_lambda_d_lambda_dif_indx(lambdas[i], lambdas[j],
					      mixing[i], mean,
					      N, values_size);
	MN_deriv[indx] += 
	  d_MN_d_real_mixing(current_lambdas[j], N)
	  *d_real_mix_d_lambda_dif_indx(lambdas[j], mixing[j], 
					lambdas[i], mixing[i], mean);
      }
    }
  }
}


void
calculate_current_lambdas(const vector<double> &lambdas,
			  const vector<double> &real_mixing,
			  const vector<double> &mixing,
			  const size_t sample_size,
			  const size_t current_reads,
			  vector<double> &current_lambdas){
  for(size_t i = 0; i < lambdas.size(); i++){
    current_lambdas[i] = exp(real_mixing[i] + log(current_reads)
			     -log(mixing[i]) - log(sample_size)
			     +log(1-exp(-lambdas[i])));
  }
}

/* delta method, returns log(var)
   remember we transform the lambdas at each N 
   and the mixings are transformed*/
static double
compute_log_var_from_fit(const vector< vector<double> > &var_matrix,
			 const vector<double> &lambdas,
			 const vector<double> &mixing,
			 const vector<double> &current_lambdas,
			 const vector<double> &real_mixing,
			 const size_t N,
			 const size_t values_size){
  double var = 0.0;
  vector<double> deriv_vec(2*lambdas.size()-1, 0.0);

  compute_deriv_MN_vec(lambdas, current_lambdas, mixing, real_mixing,
		       N, values_size, deriv_vec);

  vector<double> first_vec;
  for(size_t i = 0; i < deriv_vec.size(); i++){
    double holding_val = 0.0;
    for(size_t j = 0; j < deriv_vec.size(); j++){
      holding_val += deriv_vec[j]*var_matrix[j][i];
    } 
    first_vec.push_back(holding_val);
  }
  for(size_t i = 0; i < deriv_vec.size(); i++)
    var += deriv_vec[i]*first_vec[i];

  return(log(var));
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
/*
static double
calculate_expected_unique(const size_t sample_size,
                          const vector<double> &current_lambdas,
                          const vector<double> &mixing){

  vector<double> M_vec;
  for(size_t i = 0; i < current_lambdas.size(); i++){
    const double mean = current_lambdas[i]/(1-exp(-current_lambdas[i]));
    M_vec.push_back(mixing[i]*sample_size/mean);
  }

  return(accumulate(M_vec.begin(), M_vec.end(), 0.0));
  } */

static void
compute_log_normalizing_constant(const double lambda,
                                 double &log_normalizing_constant){
  log_normalizing_constant -= log(1-exp(-lambda)) + lambda;

}


static void
compute_normalized_lambda_vec(const double current_lambda,
			      const size_t upper_lim,
			      vector<double> &normalized_lambdas){

  normalized_lambdas[0] = current_lambda;
  for(size_t i = 1; i < upper_lim; i++)
    normalized_lambdas[i] = normalized_lambdas[i-1]*current_lambda/(i+1);

}
  
/*correct the convolution */
static void
set_initial_Qs_to_zero(const size_t iter,
                       vector<double> &current_Q){
  for(size_t i = 0; i < iter; i++)
    current_Q[i] = 0;
}

static void
set_back_Qs_to_zero(const size_t upper_lim,
                    vector<double> &current_Q){
  for(size_t i = upper_lim; i < current_Q.size(); i++)
	current_Q[i] = 0.0;
}

void
compute_current_Q_by_fft(const size_t upper_lim2,
	                 vector<double> &prev_Q,
                         vector<double> &normalized_lambdas,
                         vector<double> &current_Q){
  /*I like to use vectors, but gsl_FFt does not */
  double lambdas_array[2*upper_lim2];
  double prev_Q_array[2*upper_lim2];
  double current_Q_array[2*upper_lim2];
  for(size_t i = 0; i < upper_lim2; i++){
    REAL(lambdas_array, i) = normalized_lambdas[i];
    IMAG(lambdas_array, i) = 0.0;
    REAL(prev_Q_array, i) = prev_Q[i];
    IMAG(prev_Q_array, i) = 0.0;
    REAL(current_Q_array, i) = 0.0;
    IMAG(current_Q_array, i) = 0.0;
  }

  gsl_complex_packed_array gsl_lambdas = lambdas_array;
  gsl_complex_packed_array gsl_prev_Q = prev_Q_array;
  gsl_complex_packed_array gsl_current_Q = current_Q_array;
  /*fft by gsl */


  gsl_fft_complex_radix2_forward(gsl_lambdas, 1, upper_lim2); 
  gsl_fft_complex_radix2_forward(gsl_prev_Q, 1, upper_lim2);
  /*multiply fft(Q) and fft(lambdas) to get fft(convolution) */
  for(size_t i = 0; i < upper_lim2; i++){
    REAL(gsl_current_Q, i) = REAL(gsl_prev_Q, i)*REAL(gsl_lambdas, i) -
      IMAG(gsl_prev_Q, i)*IMAG(gsl_lambdas, i);
    IMAG(gsl_current_Q, i) = REAL(gsl_prev_Q, i)*IMAG(gsl_lambdas, i) +
      REAL(gsl_lambdas, i)*IMAG(gsl_prev_Q, i);
  }

  /*invert fft to set current_Q = convolution*/
  gsl_fft_complex_radix2_inverse(gsl_current_Q, 1, upper_lim2);
  for(size_t i = 0; i < upper_lim2; i++){
    current_Q[i] = REAL(gsl_current_Q, i);
  }
  
}

/*treat lambdas as coming from random dist as approx*/  
static double
get_rand_lambda(const vector<double> &lambdas,
		const vector<double> &mixing, 
		const gsl_rng *rng){
  double test_val  = 0.0;
  double return_lambda;
  const double unif = gsl_rng_uniform(rng);
  for(size_t i = 0; i < mixing.size(); i++){
    test_val += mixing[i];
    if(test_val > unif){
      return_lambda = lambdas[i];
      break;
    }
  }
  return(return_lambda);
}



static double
compute_log_var_from_sampling(const vector<double> &lambdas,
			  const vector<double> &current_mixing,
			  const vector<double> &real_mixing,
			  const size_t N,
			  const double tol){
  vector<double> log_expected_MN;
  for(size_t i = 0; i < lambdas.size(); i++){
    log_expected_MN.push_back(log(N) + log(real_mixing[i])
			      +log(1-exp(-lambdas[i])) 
			      -log(lambdas[i]));
  }
  const double expected_MN = exp(log_sum_log_vec(log_expected_MN,
				    log_expected_MN.size()));

  size_t upper_lim2 = 1;
  while((static_cast<double>(upper_lim2)/static_cast<double>(N))
        < 1){
    upper_lim2 = upper_lim2*2;
  }
  vector<double> prev_Q(upper_lim2, 0.0);
  prev_Q[0] = 1;

  vector<double> inverse_cdf_sum_ztp;
  
  double log_normalizing_constant = 0.0;
  size_t column_count = 0;
  vector<double> current_lambda_vec; 
 
  double test_val = std::numeric_limits<double>::max();
  size_t indx = 0;

  vector<double> log_var_vec;

  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  srand(time(0) + getpid());
  gsl_rng_set(rng, rand());

  while(test_val > tol){
    vector<double> normalized_lambda_vec(upper_lim2, 0.0);
    vector<double> current_Q(upper_lim2, 0.0);
    const double current_lambda = get_rand_lambda(lambdas,
						  current_mixing, rng);
  
    compute_normalized_lambda_vec(current_lambda, 
				  N, normalized_lambda_vec);
    compute_current_Q_by_fft(upper_lim2, prev_Q,
                               normalized_lambda_vec, current_Q);
    set_initial_Qs_to_zero(column_count, current_Q);
    set_back_Qs_to_zero(N, current_Q);
    compute_log_normalizing_constant(current_lambda, 
				     log_normalizing_constant);
    inverse_cdf_sum_ztp.push_back(accumulate(current_Q.begin(),
					     current_Q.begin()+N-1,
					     0.0)
				  *exp(log_normalizing_constant));
    test_val = inverse_cdf_sum_ztp.back();
    
    if(!finite(test_val)){
      test_val = 0.0;
    }

    if(finite(inverse_cdf_sum_ztp.back()) && indx > 1){
      double inverse_prob = 
	inverse_cdf_sum_ztp[indx-2]
	-inverse_cdf_sum_ztp[indx-1];
      if(inverse_prob < 0){
	inverse_prob = fabs(inverse_prob);
      }
      double log_var_val = 2*log(fabs(indx-2-expected_MN)) 
	+ log(inverse_prob);
      if(finite(log_var_val)){
	log_var_vec.push_back(log_var_val);
      }
    }
    indx++;

    prev_Q.swap(current_Q);
  }

  return(log_sum_log_vec(log_var_vec, log_var_vec.size()));

}

static void
calculate_current_mixings(const vector<double> &real_mixing,
                          const vector<double> &current_lambdas,
			  const size_t N,
                          vector<double> &current_mixing){

  vector<double> log_denom_vec;
  for(size_t i = 0; i < current_lambdas.size(); i++){
    current_mixing[i] = log(1-exp(-current_lambdas[i])) + log(N) 
      + log(real_mixing[i]) - log(current_lambdas[i]);
    log_denom_vec.push_back(current_mixing[i]);
  }
  const double log_denom = log_sum_log_vec(log_denom_vec, 
                                     log_denom_vec.size());

  for(size_t i = 0; i < current_mixing.size(); i++)
    current_mixing[i] = exp(current_mixing[i] - log_denom);

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
    OptionParser opt_parse("poisson_estimation_extrapolation", 
                           "fit Poisson mixture and extrapolate complexity",
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
    vector<ZTP_mixture> mixture_vec;

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
      score_vec.push_back(score);
      mixture_vec.push_back(ZTP_mix);
			   
    }

    size_t opt_num = select_number_mixtures(score_vec);
    size_t opt_num_mixs = step_size*(opt_num) + min_number_mixtures;

    ZTP_mixture opt_mix = mixture_vec[opt_num];
    vector<ZTP> opt_distros = opt_mix.get_distros();
    vector<double> opt_lambdas;
    for(size_t i = 0; i < opt_distros.size(); i++)
      opt_lambdas.push_back(opt_distros[i].get_lambda());

    vector<double> opt_mixing = opt_mix.get_mixing();
    vector< vector<double> > Fisher_info = opt_mix.get_Fish_info();

    vector<double> expected_MN;
    for(size_t i = step_btwn_extra; i <= extrapolation_size;
	i += step_btwn_extra){
      expected_MN.push_back(opt_mix.expected_inverse_sum(values_size, i));
    }    


    if(bootstraps == 0){
      invert_Fisher_info(Fisher_info);


      vector<double> current_lambdas(opt_lambdas);
      vector<double> real_mixing(opt_mixing);
      calculate_real_mixing(opt_lambdas, opt_mixing, real_mixing);

      if(sampling_error == false){
	ostream* out = (outfile.empty()) ? 
	  &std::cout : new std::ofstream(outfile.c_str());
	size_t time_step = step_btwn_extra;
	for(size_t i = 0; i < expected_MN.size(); i++){
	  calculate_current_lambdas(opt_lambdas, real_mixing, opt_mixing, 
				    values_size, time_step, current_lambdas);


	  double log_var_fit = 
	    compute_log_var_from_fit(Fisher_info, opt_lambdas,opt_mixing,
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
	vector<double> current_mixing(opt_mixing);
	double var_samp = 0.0;

	for(size_t i = 0; i < expected_MN.size(); i++){
	  calculate_current_lambdas(opt_lambdas, real_mixing, opt_mixing, 
				    values_size,time_step, current_lambdas);

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
	    exp(compute_log_var_from_fit(Fisher_info, opt_lambdas,opt_mixing,
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
      double holding_score =0.0;
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

	holding_score = ZTP_mix.EM_mix_resolve(sample_hist, tolerance, 
					       max_iter);
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
      size_t time_step = step_btwn_extra;
      for(size_t i = 0; i < bootstrap_MN[0].size(); i++){
	*out << time_step << "\t";
	for(size_t j = 0; j < bootstrap_MN.size(); j++){
	  *out << bootstrap_MN[j][i] << "\t";
	}
	*out << endl;
	time_step += step_btwn_extra;

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
