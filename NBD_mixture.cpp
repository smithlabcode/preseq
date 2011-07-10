/*    NBD_mixture:
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

#include "rmap_utils.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

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
using std::cin;

double 
log_sum_log_vec(const vector<double> &vals, size_t limit){
  const size_t max_idx = 
    max_element(vals.begin(), vals.begin() + limit) - 
    vals.begin();
  const double max_val = vals[max_idx];
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(finite(sum)); 
// abort if the sum is infinte //
#endif
    }
  }
  return(max_val + log(sum));
}

const double NBD::max_allowed_alpha = 1000;
const double NBD::min_allowed_alpha = 1e-20;
const double NBD::tolerance = 1e-10; 

void
NBD::set_helpers() {
  n_helper = 1/alpha;
  p_helper = n_helper/(n_helper + mu);
  n_log_p_minus_lngamma_n_helper = 
    n_helper*log(p_helper) - gsl_sf_lngamma(n_helper);
  log_q_helper = log(1 - p_helper);
  //TODO: check these!!!!!!!!!!!!!!!!  they are correct.
}

double 
NBD::operator()(const int val) const {
  const double P = (gsl_sf_lngamma(val + n_helper) - 
		    gsl_sf_lnfact(static_cast<size_t>(val))) + 
    n_log_p_minus_lngamma_n_helper + val*log_q_helper;
  if (!finite(P)) return -40;
  return P;
}




double
NBD::score_fun_first_term(const vector<size_t> &vals_hist,
			  const double a_mid){
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
	inner_sum += j/(1 + a_mid*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum;
}

double
NBD::alpha_score_function(const vector<size_t> &vals_hist,
			  const double mean,
			  const double a_mid,
			  const double vals_count){
  const double one_plus_alpha_mu = 1 + a_mid*mean;
  return (score_fun_first_term(vals_hist, a_mid)/vals_count + 
	  (log(one_plus_alpha_mu)/a_mid - mean)/a_mid); 
}

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

static inline double
compute_mean(const vector<size_t> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(),
				   vals_hist.end(), 0));
  double mean = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}

void 
NBD::estim_params(const vector<size_t> &vals_hist){
  mu = compute_mean(vals_hist);
  //mu= (1/n)sum(x), accumulate takes the sum of vals.begin
  
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), 
				   vals_hist.end(), 0));
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  //bisection, but what happened to the terms involving the gamma func? See Zhang et al. top of page 7
  set_helpers();
}

double
NBD::score_fun_first_term(const vector<double> &vals_hist,
			  const double a_mid) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
	inner_sum += j/(1 + a_mid*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum; 
}

double
NBD::alpha_score_function(const vector<double> &vals_hist, 
			  const double mean,
			  const double a_mid,
			  const double vals_count){
  const double one_plus_alpha_mu = 1 + a_mid*mean;
  return (score_fun_first_term(vals_hist, a_mid)/vals_count + 
	  (log(one_plus_alpha_mu)/a_mid - mean)/a_mid); 
}


static inline double
compute_mean(const vector<double> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(),
				   vals_hist.end(), 0.0));
  double mean = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}
    

void 
NBD::estim_params(const vector<size_t> &vals_hist,
		  const vector<double> &probs){
  vector<double> pseudo_hist(vals_hist.size(), 0.0);
  for(size_t i = 0; i < vals_hist.size(); i++){
    pseudo_hist[i] = vals_hist[i]*probs[i];
  }
  mu = compute_mean(pseudo_hist);
  
  const double pseudo_size = 
    accumulate(pseudo_hist.begin(), pseudo_hist.end(), 0.0);
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(pseudo_hist,mu, a_mid, pseudo_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  //bisection, but what happened to the terms involving the gamma func? See Zhang et al. top of page 7
  set_helpers();
}

double
NBD::log_pdf(const size_t val){
  double holding_val = 0.0;
  if(val > 0){
    for(size_t j =0; j < val; ++j){
      holding_val += log(1+alpha*j);
    }
    holding_val -= gsl_sf_lngamma(val+1);
  }
  return(holding_val + val*log(mu) - (val+1/alpha)*log(1+alpha*mu));
}

double
NBD::log_L(const vector<size_t> &vals_hist){
  double log_L = 0;;
  for(size_t i = 0; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += log(1+alpha*j);
    if(i > 0)
      holding_val -= gsl_sf_lngamma(i+1);
    log_L += vals_hist[i]*(holding_val + i*log(mu) 
			   -(i + 1/alpha)*log(1 + alpha*mu));
  }
    
  return(log_L);
}

const double ZTNBD::max_allowed_alpha = 1000;
const double ZTNBD::min_allowed_alpha = 1e-20;
const double ZTNBD::tolerance = 1e-10; 

double
ZTNBD::expected_zeros(const double pseudo_size){
  const double alpha = get_alpha();
  const double mu = get_mu();
  const double prob_zero = pow(1+alpha*mu, -1/alpha);
  const double expected_zeros = pseudo_size*(prob_zero/(1-prob_zero));
  return(expected_zeros);
}

double 
ZTNBD::trunc_log_pdf(const size_t val){
  double return_val = 0.0;
  double holding_val = 0.0;
  const double prob_zero = exp(-log(1+alpha*mu)/alpha);
  if(val > 0){
    for(size_t j =0; j < val; j++){
      holding_val += log(1+alpha*j);
    }
    holding_val -= gsl_sf_lngamma(val+1);
  }
  if(val > 0)
    return_val  = holding_val +val*log(mu) - (val+1/alpha)*log(1+alpha*mu)
    -log(1-prob_zero);
  else
    return_val  = -log(1+alpha*mu)/alpha;
  
  return(return_val);
}

double
ZTNBD::log_pdf(const size_t val){
  double holding_val = 0.0;
  if(val > 0){
    for(size_t j =0; j < val; j++){
      holding_val += log(1+alpha*j);
    }
    holding_val -= gsl_sf_lngamma(val+1);
  }
  return(holding_val + val*log(mu) - (val+1/alpha)*log(1+alpha*mu));
}

void
ZTNBD::set_helpers() {
  n_helper = 1/alpha;
  p_helper = n_helper/(n_helper + mu);
  n_log_p_minus_lngamma_n_helper = 
    n_helper*log(p_helper) - gsl_sf_lngamma(n_helper);
  log_q_helper = log(1 - p_helper);
  //TODO: check these!!!!!!!!!!!!!!!!  they are correct.
}

double 
ZTNBD::operator()(const int val) const {
  const double P = (gsl_sf_lngamma(val + n_helper) - 
		    gsl_sf_lnfact(static_cast<size_t>(val))) + 
    n_log_p_minus_lngamma_n_helper + val*log_q_helper;
  if (!finite(P)) return -40;
  return P;
}




double
ZTNBD::score_fun_first_term(const vector<size_t> &vals_hist,
			    const double a_mid){
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
	inner_sum += j/(1 + a_mid*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum;
}

double
ZTNBD::alpha_score_function(const vector<size_t> &vals_hist,
			    const double mean,
			    const double a_mid,
			    const double vals_count){
  const double one_plus_alpha_mu = 1 + a_mid*mean;
  return (score_fun_first_term(vals_hist, a_mid)/vals_count + 
	  (log(one_plus_alpha_mu)/a_mid - mean)/a_mid); 
}


void 
ZTNBD::estim_params(const vector<size_t> &vals_hist){
  mu = compute_mean(vals_hist);
  //mu= (1/n)sum(x), accumulate takes the sum of vals.begin
  
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), 
				   vals_hist.end(), 0));
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  //bisection, but what happened to the terms involving the gamma func? See Zhang et al. top of page 7
  set_helpers();
}

double
ZTNBD::score_fun_first_term(const vector<double> &vals_hist,
			  const double a_mid) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
	inner_sum += j/(1 + a_mid*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum; 
}

double
ZTNBD::alpha_score_function(const vector<double> &vals_hist, 
			  const double mean,
			  const double a_mid,
			  const double vals_count){
  const double one_plus_alpha_mu = 1 + a_mid*mean;
  return (score_fun_first_term(vals_hist, a_mid)/vals_count + 
	  (log(one_plus_alpha_mu)/a_mid - mean)/a_mid); 
}


void 
ZTNBD::estim_params(const vector<size_t> &vals_hist,
                    const vector<double> &probs){
  vector<double> pseudo_hist(vals_hist.size(), 0.0);
  for(size_t i = 0; i < vals_hist.size(); i++){
    pseudo_hist[i] = vals_hist[i]*probs[i];
  }
  mu = compute_mean(pseudo_hist);
  
  const double pseudo_size = 
    accumulate(pseudo_hist.begin(), pseudo_hist.end(), 0.0);
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(pseudo_hist,mu, a_mid, pseudo_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  //bisection, but what happened to the terms involving the gamma func? See Zhang et al. top of page 7
  set_helpers();
}


double
ZTNBD::trunc_log_L(const vector<size_t> &vals_hist){
  double log_L = 0.0;
  const double prob_zero = exp(-log(1+alpha*mu)/alpha);
  for(size_t i = 1; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += log(1+alpha*j);
 
    log_L += vals_hist[i]*(holding_val-gsl_sf_lngamma(i+1)
			   +i*log(mu)-(i+1/alpha)*log(1+alpha*mu)
			   -log(1-prob_zero));
  }
    
  return(log_L);
}

double 
ZTNBD::trunc_pval(const size_t val){
  double pval = 1.0;
  for(size_t i = 1; i < val; i++){  
    pval -= exp(trunc_log_pdf(i));
  }
  return(pval);
}

double
ZTNBD::EM_estim_params(const double tol, const size_t max_iter,
		      vector<size_t> &vals_hist){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(),0));
  double error = numeric_limits<double>::max();
  double prev_score = numeric_limits<double>::max();
  double score = 0.0;
  for(size_t i = 0; i < max_iter; i++){
    vals_hist[0] = expected_zeros(vals_size);
    estim_params(vals_hist);
    vals_hist[0] = 0;
    score = trunc_log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol)
      break;
    prev_score = score;
  }
  return(trunc_log_L(vals_hist));
}

double 
ZTNBD::expected_inverse_sum(const double mean,
                            const size_t sample_size,
                            const size_t sum){
  const double alph = get_alpha();
  const double m = get_mu();
  const double prob_zero = exp(-log(1+alph*m)/alph);
  return(exp(log(sample_size) - log(1-prob_zero) + 
	     log(1-exp(-log(1+exp(log(alph)+log(sum)-log(sample_size)
				  +log(m) - log(mean)))/alph))));
}

double 
NBD_mixture::log_L(const vector<size_t> &vals_hist){
  double logL = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    double inner_sum = 0.0;
    vector<double> workspace_vec(distros.size(), 0.0);
    for(size_t j = 0; j < distros.size(); j++)
      workspace_vec[j] = log(mixing[j]) + 
                         distros[j].log_pdf(i);
    
    inner_sum = log_sum_log_vec(workspace_vec, workspace_vec.size());
    logL += vals_hist[i]*inner_sum;
  }
  return(logL);
}

static inline double
d_log_L_d_mu(const double mu,
             const double alpha,
             const size_t value){
  const double val = static_cast<double>(value);
  return(val/mu - (1+alpha*val)/(1+alpha*mu));
}

static inline double
d_log_L_d_alpha(const double mu,
                const double alpha,
                const size_t value){
  const double val = static_cast<double>(value);
  double inner_sum = 0.0;
  for(size_t i = 0; i < value; i++)
    inner_sum += i/(1+alpha*i);
  return(inner_sum + log(1+alpha*mu)/(alpha*alpha)
         -mu*(1+alpha*val)/(alpha+alpha*alpha*mu));
}

static inline double
d2_log_L_d_mu2(const double mu,
               const double alpha,
               const size_t value){
  const double val = static_cast<double>(value);
  return(-val/(mu*mu) + alpha*(1+alpha*val)/((1+alpha*mu)*(1+alpha*mu)));
}

static inline double
d2_log_L_mixed(const double mu,
               const double alpha,
               const size_t value){
  const double val = static_cast<double>(value);
  return((mu-val)/((1+alpha*mu)*(1+alpha*mu)));
}

static inline double
d2_log_L_d_alpha2(const double mu,
                  const double alpha,
                  const size_t value){
  const double val = static_cast<double>(value);
  double inner_sum = 0.0;
  for(size_t i = 0; i < value; i++)
    inner_sum += i*i/((1+alpha*i)*(1+alpha*i));
  return(-inner_sum - 2*log(1+alpha*mu)/(alpha*alpha*alpha)
         +(2*mu+3*alpha*mu*mu
           +alpha*alpha*mu*mu*val)/(alpha*alpha*(1+alpha*mu)*(1+alpha*mu)));
}

// For how to compute Fisher info, see Louis(1982)

double
NBD_obs_fish_info_2nd_theta(const vector< vector<double> > &probs,
                            const vector<size_t> &vals_hist,
                            const vector<double> &mixing,
                            const size_t indx){
  const double theta = mixing[indx];
  const double last_theta = mixing.back();
  vector<double> terms_vec;

  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double prob = probs[indx][i];
      const double last_prob = probs[probs.size()-1][i];
      //-E(2nd deriv)
      terms_vec.push_back(hist_val*(prob/(theta*theta)
                                    +last_prob/(last_theta*last_theta)));
      //-E(1st deriv^2)
      terms_vec.push_back(-hist_val*(prob/(theta*theta)
                                     +last_prob/(last_theta*last_theta)));
      //E(1st deriv)^2
      const double first_deriv = prob/theta - last_prob/last_theta;
      terms_vec.push_back(hist_val*first_deriv*first_deriv);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double 
NBD_obs_fish_info_mixed_theta(const vector< vector<double> > &probs,
			      const vector<size_t> &vals_hist,
			      const vector<double> &mixing,
			      const size_t indx1,
			      const size_t indx2){
  vector<double> terms_vec;
  const double theta1 = mixing[indx1];
  const double theta2 = mixing[indx2];
  const double last_theta = mixing.back();
  for(size_t i =0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double prob1 = probs[indx1][i];
      const double prob2 = probs[indx2][i];
      const double last_prob = probs[probs.size()-1][i];
      //-E(2nd deriv)
      terms_vec.push_back(hist_val*last_prob/(last_theta*last_theta));
      //-E(1stderiv_1*1st_deriv2)
      terms_vec.push_back(-hist_val*last_prob/(last_theta*last_theta));
      //E(1stderiv_1)*E(1st_deriv_2)
      const double first_deriv1 = prob1/theta1-last_prob/last_theta;
      const double first_deriv2 = prob2/theta2-last_prob/last_theta;
      terms_vec.push_back(hist_val*first_deriv1*first_deriv2);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}



double
NBD_obs_fish_info_2nd_mu(const vector< vector<double> > &probs,
			 const vector<size_t> &vals_hist,
			 const vector<NBD> &distros,
			 const size_t indx){
  vector<double> terms_vec;
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();

  for(size_t i =0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv)
      const double second_deriv = d2_log_L_d_mu2(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*second_deriv);
      //-E(1st_deriv^2)
      const double first_deriv = d_log_L_d_mu(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*first_deriv*first_deriv);
      //E(1st_deriv)^2
      terms_vec.push_back(hist_val*prob*prob*first_deriv*first_deriv);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double 
NBD_obs_fish_info_mixed_mu(const vector< vector<double> > &probs,
			   const vector<size_t> &vals_hist,
			   const vector<NBD> &distros,
			   const size_t indx1,
			   const size_t indx2){
  vector<double> terms_vec;
  const double mu1 = distros[indx1].get_mu();
  const double mu2 = distros[indx2].get_mu();
  const double alpha1 = distros[indx1].get_alpha();
  const double alpha2 = distros[indx2].get_alpha();
  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob1 = probs[indx1][i];
      const double prob2 = probs[indx2][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv) = 0, E(1st_deriv_1*1st_deriv_2) = 0
      //E(1st_deriv_1)*E(1st_deriv_2)
      const double first_deriv1 = d_log_L_d_mu(mu1,alpha1,i);
      const double first_deriv2 = d_log_L_d_mu(mu2, alpha2, i);
      terms_vec.push_back(hist_val*prob1*prob2*first_deriv1*first_deriv2);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
NBD_obs_fish_info_2nd_alpha(const vector< vector<double> > &probs,
			    const vector<size_t> &vals_hist,
			    const vector<NBD> &distros,
			    const size_t indx){
  vector<double> terms_vec;

  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();

  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv)
      const double second_deriv = d2_log_L_d_alpha2(mu, alpha, i);
      const double first_deriv = d_log_L_d_alpha(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*second_deriv);
      //-E(1st_deriv^2)
      terms_vec.push_back(-hist_val*prob*first_deriv*first_deriv);
      //E(1st_deriv)^2
      terms_vec.push_back(hist_val*prob*prob*first_deriv*first_deriv);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
NBD_obs_fish_info_mixed_alpha(const vector< vector<double> > &probs,
			      const vector<size_t> &vals_hist,
			      const vector<NBD> &distros,
			      const double indx1,
			      const double indx2){
  vector<double> terms_vec;
  const double alpha1 = distros[indx1].get_alpha();
  const double alpha2 = distros[indx2].get_alpha();
  const double mu1 = distros[indx1].get_mu();
  const double mu2 = distros[indx2].get_mu();
  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob1 = probs[indx1][i];
      const double prob2 = probs[indx2][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0, E(1st_deriv_1*1st_deriv_2) = 0
      const double first_deriv1 = d_log_L_d_alpha(mu1, alpha1, i);
      const double first_deriv2 = d_log_L_d_alpha(mu2, alpha2, i);
      //E(1st_deriv1)*E(1st_deriv2)
      terms_vec.push_back(hist_val*prob1*prob2*first_deriv1*first_deriv2);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}


double
NBD_obs_fish_info_dtheta_dmu_same_indx(const vector< vector<double> > &probs,
				       const vector<size_t> &vals_hist,
				       const vector<NBD> &distros,
				       const vector<double> &mixing,
				       const size_t indx){
  const double last_theta = mixing.back();
  const double theta = mixing[indx];
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();

  vector<double> terms_vec;
  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double last_prob = probs[probs.size()-1][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      // E(2nd_deriv) = 0
      //-E(1st_deriv_theta*1st_deriv_mu)
      const double first_deriv = d_log_L_d_mu(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*first_deriv*(1/theta));
      // E(1st_deriv_theta)*E(1st_deriv_mu)
      terms_vec.push_back(hist_val*prob*first_deriv*(prob/theta
						     -last_prob/last_theta));
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
NBD_obs_fish_info_dtheta_dmu_diff_indx(const vector< vector<double> > &probs,
				       const vector<size_t> &vals_hist,
				       const vector<NBD> &distros,
				       const vector<double> &mixing,
				       const size_t indx_mix,
				       const size_t indx_mu){
  const double theta = mixing[indx_mix];
  const double last_theta = mixing.back();
  const double mu = distros[indx_mu].get_mu();
  const double alpha = distros[indx_mu].get_alpha();

  vector<double> terms_vec;

  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob_theta = probs[indx_mix][i];
      const double last_prob = probs[probs.size()-1][i];
      const double prob_mu = probs[indx_mu][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0, E(1st_deriv_theta*1st_deriv_mu) = 0
      const double first_deriv_mu = d_log_L_d_mu(mu, alpha, i);
      const double E_first_deriv_theta = 
	prob_theta/theta - last_prob/last_theta;
      //E(1st_deriv_theta*1st_deriv_mu)
      terms_vec.push_back(hist_val*prob_mu*first_deriv_mu*E_first_deriv_theta);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
NBD_obs_fish_info_dtheta_d_last_mu(const vector< vector<double> > &probs,
				   const vector<size_t> &vals_hist,
				   const vector<NBD> &distros,
				   const vector<double> &mixing,
				   const size_t indx){
  const double mu = distros[distros.size()-1].get_mu();
  const double alpha = distros[distros.size()-1].get_alpha();
  const double last_theta = mixing.back();
  const double theta = mixing[indx];

  vector<double> terms_vec;

  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double last_prob = probs[probs.size()-1][i];
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0
      const double first_deriv_mu = d_log_L_d_mu(mu, alpha, i);
      //-E(1st_deriv_theta*1st_deriv_mu)
      terms_vec.push_back(hist_val*last_prob*first_deriv_mu/last_theta);
      //E(1st_deriv_theta)*E(1st_deriv_mu)
      const double first_deriv_theta = prob/theta - last_prob/last_theta;
      terms_vec.push_back(hist_val*last_prob*first_deriv_mu*first_deriv_theta);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}
      
double
NBD_obs_fish_info_dtheta_dalpha_same_indx(const vector< vector<double> > &probs,
					  const vector<size_t> &vals_hist,
					  const vector<NBD> &distros,
					  const vector<double> &mixing,
					  const size_t indx){
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();
  const double theta = mixing[indx];
  const double last_theta = mixing.back();

  vector<double> terms_vec;
  for(size_t i  = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double last_prob = probs[probs.size()-1][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      //-E(1st_deriv_theta*1st_deriv_alpha)
      terms_vec.push_back(-hist_val*prob*first_deriv_alpha/theta);
      const double first_deriv_theta = prob/theta - last_prob/last_theta;
      //E(1st_deriv_theta)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob*first_deriv_alpha*first_deriv_theta);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
NBD_obs_fish_info_dtheta_dalpha_diff_indx(const vector< vector<double> > &probs,
					  const vector<size_t> &vals_hist,
					  const vector<NBD> &distros,
					  const vector<double> &mixing,
					  const size_t indx_mix,
					  const size_t indx_alpha){
  const double alpha = distros[indx_alpha].get_alpha();
  const double mu = distros[indx_alpha].get_mu();
  const double theta = mixing[indx_mix];
  const double last_theta = mixing.back();

  vector<double> terms_vec;

  for(size_t  i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob_alpha = probs[indx_alpha][i];
      const double prob_theta = probs[indx_mix][i];
      const double last_prob = probs[probs.size()-1][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      //E(2nd_deriv) = 0, E(1st_deriv_theta*1st_deriv_alpha) = 0
      const double first_deriv_theta = prob_theta/theta - last_prob/last_theta;
      //E(1st_deriv_theta)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob_alpha*first_deriv_alpha*first_deriv_theta);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double 
NBD_obs_fish_info_dtheta_d_last_alpha(const vector< vector<double> > &probs,
				      const vector<size_t> &vals_hist,
				      const vector<NBD> &distros,
				      const vector<double> &mixing,
				      const size_t indx){
  const double alpha = distros[distros.size()-1].get_alpha();
  const double mu = distros[distros.size()-1].get_mu();
  const double last_theta = mixing.back();
  const double theta = mixing[indx];
  vector<double> terms_vec;

  for(size_t  i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double last_prob = probs[probs.size()-1][i];
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      //-E(1st_deriv_theta*1st_deriv_alpha)
      terms_vec.push_back(hist_val*last_prob*first_deriv_alpha/last_theta);
      const double first_deriv_theta = prob/theta - last_prob/last_theta;
      //E(1st_deriv_theta)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*last_prob*first_deriv_alpha*first_deriv_theta);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
NBD_obs_fish_info_dmu_dalpha_same_indx(const vector< vector<double> > &probs,
				       const vector<size_t> &vals_hist,
				       const vector<NBD> &distros,
				       const size_t indx){
  const double alpha = distros[indx].get_alpha();
  const double mu = distros[indx].get_mu();

  vector<double> terms_vec;

  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double second_mix_deriv = d2_log_L_mixed(mu, alpha, i);
      const double first_deriv_mu = d_log_L_d_mu(mu, alpha, i);
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv)
      terms_vec.push_back(-hist_val*prob*second_mix_deriv);
      //-E(1st_deriv_mu*1st_deriv_alpha)
      terms_vec.push_back(-hist_val*prob*first_deriv_mu*first_deriv_alpha);
      //E(1st_deriv_mu)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob*prob*first_deriv_mu*first_deriv_alpha);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
NBD_obs_fish_info_dmu_dalpha_diff_indx(const vector< vector<double> > &probs,
				       const vector<size_t> &vals_hist,
				       const vector<NBD> &distros,
				       const size_t indx_mu,
				       const size_t indx_alpha){
  const double alpha1 = distros[indx_mu].get_alpha();
  const double mu1 = distros[indx_mu].get_mu();
  const double alpha2 = distros[indx_alpha].get_alpha();
  const double mu2 = distros[indx_alpha].get_mu();


  vector<double> terms_vec;
  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob_mu = probs[indx_mu][i];
      const double prob_alpha = probs[indx_alpha][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double deriv_mu = d_log_L_d_mu(mu1, alpha1, i);
      const double deriv_alpha = d_log_L_d_alpha(mu2, alpha2, i);
      //E(2nd_deriv) = 0, E(1st_deriv_mu*1st_deriv_alpha) = 0
      //E(1st_deriv_mu)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob_mu*prob_alpha*deriv_mu*deriv_alpha);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}






void
NBD_mixture::compute_Fisher_info(const vector< vector<double> > &probs,
				 const vector<size_t> &vals_hist,
				 const vector<double> &mixing){
  const vector<NBD> distros = get_distros();
  const size_t num_mixs=mixing.size();
  vector< vector<double> > holding_fish(3*num_mixs-1, 
					vector<double>(3*num_mixs-1, 0.0));
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = i; j < num_mixs-1; j++){
      if(j == i){
	holding_fish[i][j] = NBD_obs_fish_info_2nd_theta(probs, vals_hist,
							 mixing, i);	
      }
      else{
	holding_fish[i][j] = NBD_obs_fish_info_mixed_theta(probs, vals_hist,
							   mixing, i, j);
	holding_fish[j][i] = holding_fish[i][j];
      }
    }
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = i; j < num_mixs; j++){
      if(j==i){
	holding_fish[i+num_mixs-1][j+num_mixs-1] = 
	  NBD_obs_fish_info_2nd_mu(probs, vals_hist, distros, i);

      }
      else{
	holding_fish[i+num_mixs-1][j+num_mixs-1] = 
	  NBD_obs_fish_info_mixed_mu(probs, vals_hist, distros, i, j);
	holding_fish[j+num_mixs-1][i+num_mixs-1] = 
	  holding_fish[i+num_mixs-1][j+num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = i; j < num_mixs; j++){
      if(j==i){
	holding_fish[i+2*num_mixs-1][j+2*num_mixs-1] = 
	  NBD_obs_fish_info_2nd_alpha(probs, vals_hist, distros, i);
      }
      else{
	holding_fish[i+2*num_mixs-1][j+2*num_mixs-1] = 
	  NBD_obs_fish_info_mixed_alpha(probs, vals_hist, distros, i, j);
	holding_fish[j+2*num_mixs-1][i+2*num_mixs-1] = 
	  holding_fish[i+2*num_mixs-1][j+2*num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = 0; j < num_mixs-1; j++){
      if(j==i){
	holding_fish[i][j+num_mixs-1] = 
	  NBD_obs_fish_info_dtheta_dmu_same_indx(probs, vals_hist,
						 distros, mixing, i);
	holding_fish[j+num_mixs-1][i] = 
	  holding_fish[i][j+num_mixs-1];
      }
      else{
	holding_fish[i][j+num_mixs-1] = 
	  NBD_obs_fish_info_dtheta_dmu_diff_indx(probs, vals_hist,
						 distros, mixing, i, j);
	holding_fish[j+num_mixs-1][i] = 
	  holding_fish[i][j+num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    holding_fish[i][2*num_mixs-2] = 
      NBD_obs_fish_info_dtheta_d_last_mu(probs, vals_hist, distros,
					 mixing, i);
    holding_fish[2*num_mixs-2][i] = 
      holding_fish[i][2*num_mixs-2];
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = 0; j < num_mixs-1; j++){
      if(j==i){
	holding_fish[i][j+2*num_mixs-1] = 
	  NBD_obs_fish_info_dtheta_dalpha_same_indx(probs, vals_hist,
						    distros, mixing, i);
	holding_fish[j+2*num_mixs-1][i] = 
	  holding_fish[i][j+2*num_mixs-1];
      }
      else{
	holding_fish[i][j+2*num_mixs-1] = 
	  NBD_obs_fish_info_dtheta_dalpha_diff_indx(probs, vals_hist,
						    distros, mixing, i, j);
	holding_fish[j+2*num_mixs-1][i] = 
	  holding_fish[i][j+2*num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    holding_fish[i][3*num_mixs-2] = 
      NBD_obs_fish_info_dtheta_d_last_alpha(probs, vals_hist,
					    distros, mixing, i);
    holding_fish[3*num_mixs-2][i] = holding_fish[i][3*num_mixs-2];
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = 0; j < num_mixs; j++){
      if(j==i){
	holding_fish[2*num_mixs-1+i][num_mixs-1+j] = 
	  NBD_obs_fish_info_dmu_dalpha_same_indx(probs, vals_hist,
						 distros, i);
	holding_fish[num_mixs-1+j][2*num_mixs-1+i] = 
	  holding_fish[2*num_mixs-1+i][num_mixs-1+j];
      }
      else{
	holding_fish[2*num_mixs-1+i][num_mixs-1+j] =
	  NBD_obs_fish_info_dmu_dalpha_diff_indx(probs, vals_hist,
						 distros, j, i);
	holding_fish[num_mixs-1+j][2*num_mixs-1+i] = 
	  holding_fish[2*num_mixs-1+i][num_mixs-1+j];
      }
    }
  }
  set_Fish_info(holding_fish);

}

void
NBD_mixture::calculate_mixing(const vector<size_t> &vals_hist,
                 const vector< vector<double> > &probs){
  double vals_size = static_cast<double>(accumulate(vals_hist.begin(),
                                                    vals_hist.end(), 0));

  for(size_t i = 0; i < mixing.size(); i++){
    double inner_sum = 0.0;
    for(size_t j = 0; j < vals_hist.size(); j ++)
      inner_sum += vals_hist[j]*probs[i][j];
    mixing[i] = inner_sum/vals_size;
    assert(finite(mixing[i]));
  }
  double mixing_sum = accumulate(mixing.begin(), mixing.end(), 0.0);
  for(size_t i = 0; i < mixing.size(); i++){
    mixing[i] = mixing[i]/mixing_sum;
    assert(mixing[i] >= 0);
  }
}



void
NBD_mixture::expectation_step(const vector<size_t> &vals_hist,
			      vector< vector<double> > &probs){
  for(size_t i = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      vector<double> log_denom_vec;
    
      for(size_t j = 0; j < distros.size(); j++){
	log_denom_vec.push_back(log(mixing[j]) +
				distros[j].log_pdf(i));
	if(!finite(log_denom_vec[j])){
	  log_denom_vec[j] = -std::numeric_limits<double>::max()/
	    static_cast<double>(distros.size()); // arbitrary
	}
      }
    
      const double log_denom = log_sum_log_vec(log_denom_vec, 
					       log_denom_vec.size());
      assert(finite(log_denom));
  
      for(size_t j = 0; j < distros.size(); j++)
	probs[j][i] = exp(log_denom_vec[j] - log_denom);
    }
  }
}


void
NBD_mixture::maximization_step(const vector<size_t> &vals_hist,
			       const vector< vector<double> > &probs){

  for(size_t i = 0; i < distros.size(); i++){
    distros[i].estim_params(vals_hist, probs[i]);
  }
}


double 
NBD_mixture::EM_resolve_mix(const vector<size_t> &vals_hist,
			    const double &tol,const size_t max_iter){
  const size_t number_states = distros.size();
  double probs_starting_val = 1/static_cast<double>(number_states);
  vector< vector<double> > probs(number_states, 
                                 vector<double>(vals_hist.size(),
                                                probs_starting_val));

  double score = 0.0;
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iter; ++i){

    expectation_step(vals_hist, probs);

    maximization_step(vals_hist, probs);

    calculate_mixing(vals_hist, probs);

    score = log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;

  }
  compute_Fisher_info(probs, vals_hist, mixing);
  return(log_L(vals_hist));
}

static inline double
d_log_1_minus_prob_zero_d_mu(const double mu, 
			     const double alpha){
  const double prob_zero = exp(-log(1+alpha*mu)/alpha);
  return(exp(-(alpha+1)*log(1+alpha*mu)/alpha)/(1-prob_zero));
}

static inline double
d_log_1_minus_prob_zero_d_alpha(const double mu,
				const double alpha){
  const double prob_zero = exp(-log(1+alpha*mu)/alpha);
  return(prob_zero*(mu/(alpha*(1+alpha*mu)) 
		    - log(1+alpha*mu)/(alpha*alpha))/(1-prob_zero));
}

static inline double
d2_log_1_minus_prob_zero_d_mu2(const double mu,
			       const double alpha){
  const double prob0 = exp(-log(1+alpha*mu)/alpha);
  const double num = -(alpha+1)*prob0 + alpha*prob0*prob0;
  const double denom = (1+alpha*mu)*(1-prob0);
  return(num/(denom*denom));
}

static inline double
d2_log_1_minus_prob_zero_mixed(const double mu,
			       const double alpha){
  const double prob_zero = exp(-log(1+alpha*mu)/alpha);
  const double denom = (1 - prob_zero)*(1 - prob_zero);
  const double num = 
    -mu*exp(-(2*alpha+1)*log(1+alpha*mu)/alpha)
    + mu*exp(-2*(alpha+1)*log(1+alpha*mu)/alpha)
    +log(1+alpha*mu)*exp(-(alpha+1)*log(1+alpha*mu)/alpha)/(alpha*alpha)
    -mu*exp(-(2*alpha+1)*log(1+alpha*mu)/alpha)/alpha;
  return(num/denom);
}

static inline double
d2_log_1_minus_prob_zero_d_alpha2(const double mu, 
				  const double alpha){
  const double prob0 = exp(-log(1+alpha*mu)/alpha);
  const double first_term =
    prob0*(2*log(1+alpha*mu)/(alpha*alpha*alpha)
	   -mu/(alpha*alpha*(1+alpha*mu))
	   -mu*(1+2*alpha*mu)/(alpha*alpha*(1+alpha*mu)*(1+alpha*mu)));
  const double inner_deriv = 
    mu/(alpha*(1+alpha*mu)) - log(1+alpha*mu)/(alpha*alpha);
  const double second_term = -prob0*inner_deriv*inner_deriv;
  return(first_term/(1-prob0) + second_term/((1-prob0)*(1-prob0)));
}

static inline double
trunc_d_log_L_d_mu(const double mu, 
		   const double alpha,
		   const size_t value){
  return(d_log_L_d_mu(mu, alpha, value)
	 -d_log_1_minus_prob_zero_d_mu(mu, alpha));
}

static inline double 
trunc_d_log_L_d_alpha(const double mu, 
		      const double alpha,
		      const double value){
  return(d_log_L_d_alpha(mu, alpha, value)
	 -d_log_1_minus_prob_zero_d_alpha(mu, alpha));
}

static inline double 
trunc_d2_log_L_d_mu2(const double mu, 
		     const double alpha,
		     const size_t value){
  return(d2_log_L_d_mu2(mu, alpha, value)
	 -d2_log_1_minus_prob_zero_d_mu2(mu, alpha));
}

static inline double
trunc_d2_log_L_mixed(const double mu,
		     const double alpha, 
		     const size_t value){
  return(d2_log_L_mixed(mu, alpha, value)
	 -d2_log_1_minus_prob_zero_mixed(mu, alpha));
}

static inline double
trunc_d2_log_L_d_alpha2(const double mu,
			const double alpha,
			const size_t value){
  return(d2_log_L_d_alpha2(mu, alpha, value)
	 -d2_log_1_minus_prob_zero_d_alpha2(mu, alpha));
}

double
ZTNBD_obs_fish_info_2nd_theta(const vector< vector<double> > &probs,
			      const vector<size_t> &vals_hist,
			      const vector<double> &mixing,
			      const double expected_zeros,
			      const size_t indx){
  const double theta = mixing[indx];
  const double last_theta = mixing.back();
  vector<double> terms_vec;

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double prob = probs[indx][i];
      const double last_prob = probs[probs.size()-1][i];
      //-E(2nd deriv)
      terms_vec.push_back(hist_val*(prob/(theta*theta)
                                    +last_prob/(last_theta*last_theta)));
      //-E(1st deriv^2)
      terms_vec.push_back(-hist_val*(prob/(theta*theta)
                                     +last_prob/(last_theta*last_theta)));
      //E(1st deriv)^2
      const double first_deriv = prob/theta - last_prob/last_theta;
      terms_vec.push_back(hist_val*first_deriv*first_deriv);
    }
  }
  /*
  //I_zeros
  const double prob = probs[indx][0];
  const double last_prob = probs[probs.size()-1][0];
  terms_vec.push_back(expected_zeros*(prob/(theta*theta)
				      +last_prob/(last_theta*last_theta)));
  terms_vec.push_back(-expected_zeros*(prob/(theta*theta)
				       +last_prob/(last_theta*last_theta)));
  const double first_deriv = prob/theta - last_prob/last_theta;
  terms_vec.push_back(expected_zeros*first_deriv*first_deriv);		      
  */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double 
ZTNBD_obs_fish_info_mixed_theta(const vector< vector<double> > &probs,
				const vector<size_t> &vals_hist,
				const vector<double> &mixing,
				const double expected_zeros,
				const size_t indx1,
				const size_t indx2){
  vector<double> terms_vec;
  const double theta1 = mixing[indx1];
  const double theta2 = mixing[indx2];
  const double last_theta = mixing.back();
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double prob1 = probs[indx1][i];
      const double prob2 = probs[indx2][i];
      const double last_prob = probs[probs.size()-1][i];
      //-E(2nd deriv)
      terms_vec.push_back(hist_val*last_prob/(last_theta*last_theta));
      //-E(1stderiv_1*1st_deriv2)
      terms_vec.push_back(-hist_val*last_prob/(last_theta*last_theta));
      //E(1stderiv_1)*E(1st_deriv_2)
      const double first_deriv1 = prob1/theta1-last_prob/last_theta;
      const double first_deriv2 = prob2/theta2-last_prob/last_theta;
      terms_vec.push_back(hist_val*first_deriv1*first_deriv2);
    }
  }
  /* 
  //-I_zeros
  const double prob1 = probs[indx1][0];
  const double prob2 = probs[indx2][0];
  const double last_prob = probs[probs.size()-1][0];
  //-E(1stderiv_1*1st_deriv2 | x = zero)
  terms_vec.push_back(-expected_zeros*last_prob/(last_theta*last_theta));
  //E(1stderiv_1)*E(1st_deriv_2)
  const double first_deriv1 = prob1/theta1-last_prob/last_theta;
  const double first_deriv2 = prob2/theta2-last_prob/last_theta;
  terms_vec.push_back(expected_zeros*first_deriv1*first_deriv2);
  */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}



double
ZTNBD_obs_fish_info_2nd_mu(const vector< vector<double> > &probs,
			   const vector<size_t> &vals_hist,
			   const vector<ZTNBD> &distros,
			   const double expected_zeros,
			   const size_t indx){
  vector<double> terms_vec;
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv)
      const double second_deriv = d2_log_L_d_mu2(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*second_deriv);
      //-E(1st_deriv^2)
      const double first_deriv = d_log_L_d_mu(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*first_deriv*first_deriv);
      //E(1st_deriv)^2
      terms_vec.push_back(hist_val*prob*prob*first_deriv*first_deriv);
    }
  }
  
  // -I_zeros = -E(1st_deriv^2 | x=0) + E(1st_deriv|x=0)^2
  const double prob = probs[indx][0];
  const double second_deriv = d2_log_L_d_mu2(mu, alpha, 0);
  terms_vec.push_back(-expected_zeros*prob*second_deriv);
  const double first_deriv = d_log_L_d_mu(mu, alpha, 0);
  terms_vec.push_back(-expected_zeros*prob*first_deriv*first_deriv);
  terms_vec.push_back(expected_zeros*prob*prob*first_deriv*first_deriv);

  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double 
ZTNBD_obs_fish_info_mixed_mu(const vector< vector<double> > &probs,
			     const vector<size_t> &vals_hist,
			     const vector<ZTNBD> &distros,
			     const double expected_zeros,
			     const size_t indx1,
			     const size_t indx2){
  vector<double> terms_vec;
  const double mu1 = distros[indx1].get_mu();
  const double mu2 = distros[indx2].get_mu();
  const double alpha1 = distros[indx1].get_alpha();
  const double alpha2 = distros[indx2].get_alpha();
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob1 = probs[indx1][i];
      const double prob2 = probs[indx2][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv) = 0, E(1st_deriv_1*1st_deriv_2) = 0
      //E(1st_deriv_1)*E(1st_deriv_2)
      const double first_deriv1 = d_log_L_d_mu(mu1,alpha1,i);
      const double first_deriv2 = d_log_L_d_mu(mu2, alpha2, i);
      terms_vec.push_back(hist_val*prob1*prob2*first_deriv1*first_deriv2);
    }
  }
  
  //-I_zeros
  const double prob1=probs[indx1][0];
  const double prob2 = probs[indx2][0];
  const double first_deriv1 = d_log_L_d_mu(mu1,alpha1,0);
  const double first_deriv2 = d_log_L_d_mu(mu2, alpha2, 0);
  terms_vec.push_back(expected_zeros*prob1*prob2*first_deriv1*first_deriv2);
   
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
ZTNBD_obs_fish_info_2nd_alpha(const vector< vector<double> > &probs,
			      const vector<size_t> &vals_hist,
			      const vector<ZTNBD> &distros,
			      const double expected_zeros,
			      const size_t indx){
  vector<double> terms_vec;

  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv)
      const double second_deriv = d2_log_L_d_alpha2(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*second_deriv);
 
    }
  }
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double first_deriv = d_log_L_d_alpha(mu, alpha, i);
      //-E(1st_deriv^2)
      terms_vec.push_back(-hist_val*prob*first_deriv*first_deriv);
      //E(1st_deriv)^2
      terms_vec.push_back(hist_val*prob*prob*first_deriv*first_deriv);
    }
  }
  
  //-I_zeros
  const double prob = probs[indx][0];
  const double second_deriv = d2_log_L_d_alpha2(mu, alpha, 0);
  terms_vec.push_back(-expected_zeros*prob*second_deriv);
  const double first_deriv = d_log_L_d_alpha(mu, alpha, 0);
  terms_vec.push_back(-expected_zeros*prob*first_deriv*first_deriv);
  terms_vec.push_back(expected_zeros*prob*prob*first_deriv*first_deriv);

  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
ZTNBD_obs_fish_info_mixed_alpha(const vector< vector<double> > &probs,
				const vector<size_t> &vals_hist,
				const vector<ZTNBD> &distros,
				const double expected_zeros,
				const double indx1,
				const double indx2){
  vector<double> terms_vec;
  const double alpha1 = distros[indx1].get_alpha();
  const double alpha2 = distros[indx2].get_alpha();
  const double mu1 = distros[indx1].get_mu();
  const double mu2 = distros[indx2].get_mu();
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob1 = probs[indx1][i];
      const double prob2 = probs[indx2][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0, E(1st_deriv_1*1st_deriv_2) = 0
      const double first_deriv1 = d_log_L_d_alpha(mu1, alpha1, i);
      const double first_deriv2 = d_log_L_d_alpha(mu2, alpha2, i);
      //E(1st_deriv1)*E(1st_deriv2)
      terms_vec.push_back(hist_val*prob1*prob2*first_deriv1*first_deriv2);
    }
  }
 
  //-I_zeros
  const double prob1= probs[indx1][0];
  const double prob2 = probs[indx2][0];
  const double first_deriv1 = d_log_L_d_alpha(mu1, alpha1, 0);
  const double first_deriv2 = d_log_L_d_alpha(mu2, alpha2, 0);
  terms_vec.push_back(expected_zeros*prob1*prob2*first_deriv1*first_deriv2);
  
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}


double
ZTNBD_obs_fish_info_dtheta_dmu_same_indx(const vector< vector<double> > &probs,
					 const vector<size_t> &vals_hist,
					 const vector<ZTNBD> &distros,
					 const vector<double> &mixing,
					 const double expected_zeros,
					 const size_t indx){
  const double last_theta = mixing.back();
  const double theta = mixing[indx];
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();

  vector<double> terms_vec;
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double last_prob = probs[probs.size()-1][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      // E(2nd_deriv) = 0
      //-E(1st_deriv_theta*1st_deriv_mu)
      const double first_deriv = d_log_L_d_mu(mu, alpha, i);
      terms_vec.push_back(-hist_val*prob*first_deriv/theta);
      // E(1st_deriv_theta)*E(1st_deriv_mu)
      terms_vec.push_back(hist_val*prob*first_deriv*(prob/theta
						     -last_prob/last_theta));
    }
  }
  /*
  //-I_zeros
  const double prob = probs[indx][0];
  const double last_prob = probs[probs.size()-1][0];
  const double first_deriv  = d_log_L_d_mu(mu, alpha, 0);
  terms_vec.push_back(-expected_zeros*prob*first_deriv/theta);
  terms_vec.push_back(expected_zeros*prob*first_deriv*(prob/theta
  					       -last_prob/last_theta));
  */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
ZTNBD_obs_fish_info_dtheta_dmu_diff_indx(const vector< vector<double> > &probs,
					 const vector<size_t> &vals_hist,
					 const vector<ZTNBD> &distros,
					 const vector<double> &mixing,
					 const double expected_zeros,
					 const size_t indx_mix,
					 const size_t indx_mu){
  const double theta = mixing[indx_mix];
  const double last_theta = mixing.back();
  const double mu = distros[indx_mu].get_mu();
  const double alpha = distros[indx_mu].get_alpha();

  vector<double> terms_vec;

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob_theta = probs[indx_mix][i];
      const double last_prob = probs[probs.size()-1][i];
      const double prob_mu = probs[indx_mu][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0, E(1st_deriv_theta*1st_deriv_mu) = 0
      const double first_deriv_mu = d_log_L_d_mu(mu, alpha, i);
      const double deriv_theta = 
	prob_theta/theta - last_prob/last_theta;
      //E(1st_deriv_theta*1st_deriv_mu)
      terms_vec.push_back(hist_val*prob_mu*first_deriv_mu*deriv_theta);
    }
  }
  /*
  //-I_zeros
  const double prob_mu = probs[indx_mu][0];
  const double last_prob = probs[probs.size()-1][0];
  const double prob_theta = probs[indx_mix][0];
  const double first_deriv_mu = d_log_L_d_mu(mu, alpha, 0);
  const double deriv_theta = prob_theta/theta - last_prob/last_theta;
  terms_vec.push_back(expected_zeros*prob_mu*first_deriv_mu*deriv_theta);
  */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
ZTNBD_obs_fish_info_dtheta_d_last_mu(const vector< vector<double> > &probs,
				     const vector<size_t> &vals_hist,
				     const vector<ZTNBD> &distros,
				     const vector<double> &mixing,
				     const double expected_zeros,
				     const size_t indx){
  const double mu = distros[distros.size()-1].get_mu();
  const double alpha = distros[distros.size()-1].get_alpha();
  const double last_theta = mixing.back();
  const double theta = mixing[indx];

  vector<double> terms_vec;

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double last_prob = probs[probs.size()-1][i];
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0
      const double first_deriv_mu = d_log_L_d_mu(mu, alpha, i);
      //-E(1st_deriv_theta*1st_deriv_mu)
      terms_vec.push_back(hist_val*last_prob*first_deriv_mu/last_theta);
      //E(1st_deriv_theta)*E(1st_deriv_mu)
      const double first_deriv_theta = prob/theta - last_prob/last_theta;
      terms_vec.push_back(hist_val*last_prob*first_deriv_mu*first_deriv_theta);
    }
  }
  /*
  //-I_zeros
  const double last_prob = probs[probs.size()-1][0];
  const double prob = probs[indx][0];
  const double first_deriv_mu = d_log_L_d_mu(mu, alpha, 0);
  const double deriv_theta = prob/theta - last_prob/last_theta;
  terms_vec.push_back(expected_zeros*last_prob*first_deriv_mu/last_theta);
  terms_vec.push_back(expected_zeros*last_prob*first_deriv_mu*deriv_theta);
   */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}
      
double
ZTNBD_obs_fish_info_dtheta_dalpha_same_indx(const vector< vector<double> > &probs,
					    const vector<size_t> &vals_hist,
					    const vector<ZTNBD> &distros,
					    const vector<double> &mixing,
					    const double expected_zeros,
					    const size_t indx){
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();
  const double theta = mixing[indx];
  const double last_theta = mixing.back();

  vector<double> terms_vec;
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double last_prob = probs[probs.size()-1][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      //-E(1st_deriv_theta*1st_deriv_alpha)
      terms_vec.push_back(-hist_val*prob*first_deriv_alpha/theta);
      const double first_deriv_theta = prob/theta - last_prob/last_theta;
      //E(1st_deriv_theta)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob*first_deriv_alpha*first_deriv_theta);
    }
  }
  /*
  //-I_zeros
  const double prob = probs[indx][0];
  const double last_prob = probs[probs.size()-1][0];
  const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, 0);
  const double first_deriv_theta  = prob/theta - last_prob/last_theta;
  terms_vec.push_back(-expected_zeros*prob*first_deriv_alpha/theta);
  terms_vec.push_back(expected_zeros*prob*first_deriv_alpha*first_deriv_theta);
  */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
ZTNBD_obs_fish_info_dtheta_dalpha_diff_indx(const vector< vector<double> > &probs,
					    const vector<size_t> &vals_hist,
					    const vector<ZTNBD> &distros,
					    const vector<double> &mixing,
					    const double expected_zeros,
					    const size_t indx_mix,
					    const size_t indx_alpha){
  const double alpha = distros[indx_alpha].get_alpha();
  const double mu = distros[indx_alpha].get_mu();
  const double theta = mixing[indx_mix];
  const double last_theta = mixing.back();

  vector<double> terms_vec;

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob_alpha = probs[indx_alpha][i];
      const double prob_theta = probs[indx_mix][i];
      const double last_prob = probs[probs.size()-1][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      //E(2nd_deriv) = 0, E(1st_deriv_theta*1st_deriv_alpha) = 0
      const double deriv_theta = prob_theta/theta - last_prob/last_theta;
      //E(1st_deriv_theta)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob_alpha*first_deriv_alpha*deriv_theta);
    }
  }
  /*
  //-I_zeros
  const double prob_alpha = probs[indx_alpha][0];
  const double prob_theta = probs[indx_mix][0];
  const double last_prob = probs[probs.size()-1][0];
  const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, 0);
  const double deriv_theta = prob_theta/theta - last_prob/last_theta;
  terms_vec.push_back(expected_zeros*prob_alpha*first_deriv_alpha*deriv_theta);
  */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double 
ZTNBD_obs_fish_info_dtheta_d_last_alpha(const vector< vector<double> > &probs,
					const vector<size_t> &vals_hist,
					const vector<ZTNBD> &distros,
					const vector<double> &mixing,
					const double expected_zeros, 
					const size_t indx){
  const double alpha = distros[distros.size()-1].get_alpha();
  const double mu = distros[distros.size()-1].get_mu();
  const double last_theta = mixing.back();
  const double theta = mixing[indx];
  vector<double> terms_vec;

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double last_prob = probs[probs.size()-1][i];
      const double prob = probs[indx][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      //E(2nd_deriv) = 0
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      //-E(1st_deriv_theta*1st_deriv_alpha)
      terms_vec.push_back(hist_val*last_prob*first_deriv_alpha/last_theta);
      const double deriv_theta = prob/theta - last_prob/last_theta;
      //E(1st_deriv_theta)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*last_prob*first_deriv_alpha*deriv_theta);
    }
  }
  /*  
  //-I_zeros
  const double last_prob = probs[probs.size()-1][0];
  const double prob = probs[indx][0];
  const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, 0);
  const double deriv_theta = prob/theta - last_prob/last_theta;
  terms_vec.push_back(expected_zeros*last_prob*first_deriv_alpha/last_theta);
  terms_vec.push_back(expected_zeros*last_prob*first_deriv_alpha*deriv_theta);
  */
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
ZTNBD_obs_fish_info_dmu_dalpha_same_indx(const vector< vector<double> > &probs,
					 const vector<size_t> &vals_hist,
					 const vector<ZTNBD> &distros,
					 const double expected_zeros,
					 const size_t indx){
  const double alpha = distros[indx].get_alpha();
  const double mu = distros[indx].get_mu();

  vector<double> terms_vec;

  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob = probs[indx][i];
      const double second_mix_deriv = d2_log_L_mixed(mu, alpha, i);
      const double first_deriv_mu = d_log_L_d_mu(mu, alpha, i);
      const double first_deriv_alpha = d_log_L_d_alpha(mu, alpha, i);
      const double hist_val = static_cast<double>(vals_hist[i]);
      //-E(2nd_deriv)
      terms_vec.push_back(-hist_val*prob*second_mix_deriv);
      //-E(1st_deriv_mu*1st_deriv_alpha)
      terms_vec.push_back(-hist_val*prob*first_deriv_mu*first_deriv_alpha);
      //E(1st_deriv_mu)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob*prob*first_deriv_mu*first_deriv_alpha);
    }
  }
  
  //-I_zeros
  const double prob = probs[indx][0];
  const double second_mix_deriv = d2_log_L_mixed(mu, alpha, 0);
  terms_vec.push_back(-expected_zeros*prob*second_mix_deriv);
  const double deriv_mu = d_log_L_d_mu(mu, alpha, 0);
  const double deriv_alpha = d_log_L_d_alpha(mu, alpha, 0);
  terms_vec.push_back(-expected_zeros*prob*deriv_mu*deriv_alpha);
  terms_vec.push_back(expected_zeros*prob*prob*deriv_mu*deriv_alpha);

  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
ZTNBD_obs_fish_info_dmu_dalpha_diff_indx(const vector< vector<double> > &probs,
					 const vector<size_t> &vals_hist,
					 const vector<ZTNBD> &distros,
					 const double expected_zeros,
					 const size_t indx_mu,
					 const size_t indx_alpha){
  const double alpha1 = distros[indx_mu].get_alpha();
  const double mu1 = distros[indx_mu].get_mu();
  const double alpha2 = distros[indx_alpha].get_alpha();
  const double mu2 = distros[indx_alpha].get_mu();

  vector<double> terms_vec;
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob_mu = probs[indx_mu][i];
      const double prob_alpha = probs[indx_alpha][i];
      const double hist_val = static_cast<double>(vals_hist[i]);
      const double deriv_mu = d_log_L_d_mu(mu1, alpha1, i);
      const double deriv_alpha = d_log_L_d_alpha(mu2, alpha2, i);
      //E(2nd_deriv) = 0, E(1st_deriv_mu*1st_deriv_alpha) = 0
      //E(1st_deriv_mu)*E(1st_deriv_alpha)
      terms_vec.push_back(hist_val*prob_mu*prob_alpha*deriv_mu*deriv_alpha);
    }
  }
  
  //-I_zeros
  const double prob_mu = probs[indx_mu][0];
  const double prob_alpha = probs[indx_alpha][0];
  const double deriv_mu = d_log_L_d_mu(mu1, alpha1, 0);
  const double deriv_alpha = d_log_L_d_alpha(mu2, alpha2, 0);
  terms_vec.push_back(expected_zeros*prob_mu*prob_alpha*deriv_mu*deriv_alpha);
  
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}






void
ZTNBD_mixture::compute_Fisher_info(const vector< vector<double> > &probs,
				   const vector<size_t> &vals_hist,
				   const vector<double> &mixing,
				   const double expected_zeros){
  const vector<ZTNBD> distros = get_distros();
  const size_t num_mixs=mixing.size();
  vector< vector<double> > holding_fish(3*num_mixs-1, 
					vector<double>(3*num_mixs-1, 0.0));
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = 0; j < num_mixs-1; j++){
      if(j == i){
	holding_fish[i][j] = ZTNBD_obs_fish_info_2nd_theta(probs, vals_hist,
							   mixing, 
							   expected_zeros, i);	
      }
      else{
	holding_fish[i][j] = ZTNBD_obs_fish_info_mixed_theta(probs, vals_hist,
							     mixing, 
							     expected_zeros, 
							     i, j);
      }
    }
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = 0; j < num_mixs; j++){
      if(j==i){
	holding_fish[i+num_mixs-1][j+num_mixs-1] = 
	  ZTNBD_obs_fish_info_2nd_mu(probs, vals_hist, distros, 
				     expected_zeros, i);

      }
      else{
	holding_fish[i+num_mixs-1][j+num_mixs-1] = 
	  ZTNBD_obs_fish_info_mixed_mu(probs, vals_hist, distros, 
				       expected_zeros, i, j);
      }
    }
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = 0; j < num_mixs; j++){
      if(j==i){
	holding_fish[i+2*num_mixs-1][j+2*num_mixs-1] = 
	  ZTNBD_obs_fish_info_2nd_alpha(probs, vals_hist, distros, 
					expected_zeros, i);
      }
      else{
	holding_fish[i+2*num_mixs-1][j+2*num_mixs-1] = 
	  ZTNBD_obs_fish_info_mixed_alpha(probs, vals_hist, distros, 
					  expected_zeros, i, j);
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = 0; j < num_mixs-1; j++){
      if(j==i){
	holding_fish[i][j+num_mixs-1] = 
	  ZTNBD_obs_fish_info_dtheta_dmu_same_indx(probs, vals_hist,
						   distros, mixing, 
						   expected_zeros, i);
	holding_fish[j+num_mixs-1][i] = 
	  holding_fish[i][j+num_mixs-1];
      }
      else{
	holding_fish[i][j+num_mixs-1] = 
	  ZTNBD_obs_fish_info_dtheta_dmu_diff_indx(probs, vals_hist,
						   distros, mixing, 
						   expected_zeros, i, j);
	holding_fish[j+num_mixs-1][i] = 
	  holding_fish[i][j+num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    holding_fish[i][2*num_mixs-2] = 
      ZTNBD_obs_fish_info_dtheta_d_last_mu(probs, vals_hist, distros,
					   mixing, expected_zeros, i);
    holding_fish[2*num_mixs-2][i] = 
      holding_fish[i][2*num_mixs-2];
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = 0; j < num_mixs-1; j++){
      if(j==i){
	holding_fish[i][j+2*num_mixs-1] = 
	  ZTNBD_obs_fish_info_dtheta_dalpha_same_indx(probs, vals_hist,
						      distros, mixing, 
						      expected_zeros, i);
	holding_fish[j+2*num_mixs-1][i] = 
	  holding_fish[i][j+2*num_mixs-1];
      }
      else{
	holding_fish[i][j+2*num_mixs-1] = 
	  ZTNBD_obs_fish_info_dtheta_dalpha_diff_indx(probs, vals_hist,
						      distros, mixing, 
						      expected_zeros, i, j);
	holding_fish[j+2*num_mixs-1][i] = 
	  holding_fish[i][j+2*num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    holding_fish[i][3*num_mixs-2] = 
      ZTNBD_obs_fish_info_dtheta_d_last_alpha(probs, vals_hist,
					      distros, mixing, 
					      expected_zeros, i);
    holding_fish[3*num_mixs-2][i] = holding_fish[i][3*num_mixs-2];
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = 0; j < num_mixs; j++){
      if(j==i){
	holding_fish[2*num_mixs-1+i][num_mixs-1+j] = 
	  ZTNBD_obs_fish_info_dmu_dalpha_same_indx(probs, vals_hist,
						   distros, expected_zeros, i);
	holding_fish[num_mixs-1+j][2*num_mixs-1+i] = 
	  holding_fish[2*num_mixs-1+i][num_mixs-1+j];
      }
      else{
	holding_fish[2*num_mixs-1+i][num_mixs-1+j] =
	  ZTNBD_obs_fish_info_dmu_dalpha_diff_indx(probs, vals_hist,
						   distros, 
						   expected_zeros, j, i);
	holding_fish[num_mixs-1+j][2*num_mixs-1+i] = 
	  holding_fish[2*num_mixs-1+i][num_mixs-1+j];
      }
    }
  }
  set_Fish_info(holding_fish);

}


double 
ZTNBD_mixture::trunc_log_L(const vector<size_t> &vals_hist){
  double logL = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    double inner_sum = 0.0;
    vector<double> workspace_vec(distros.size(), 0.0);
    for(size_t j = 0; j < distros.size(); j++)
      workspace_vec[j] = log(mixing[j]) + 
                         distros[j].trunc_log_pdf(i);
    
    inner_sum = log_sum_log_vec(workspace_vec, workspace_vec.size());
    logL += vals_hist[i]*inner_sum;
  }
  return(logL);
}

void
ZTNBD_mixture::trunc_calculate_mixing(const vector<size_t> &vals_hist,
				      const vector< vector<double> > &probs){
  double vals_size = static_cast<double>(accumulate(vals_hist.begin(),
                                                    vals_hist.end(), 0));

  for(size_t i = 0; i < mixing.size(); i++){
    double inner_sum = 0.0;
    for(size_t j = 1; j < vals_hist.size(); j ++)
      inner_sum += vals_hist[j]*probs[i][j];
    mixing[i] = inner_sum/vals_size;
    assert(finite(mixing[i]));
  }
  double mixing_sum = accumulate(mixing.begin(), mixing.end(), 0.0);
  for(size_t i = 0; i < mixing.size(); i++){
    mixing[i] = mixing[i]/mixing_sum;
    assert(mixing[i] >= 0);
  }
}



void
ZTNBD_mixture::trunc_expectation_step(const vector<size_t> &vals_hist,
				      vector< vector<double> > &probs){
  for(size_t i = 0; i < vals_hist.size(); i++){
    if((i == 0) || (vals_hist[i] > 0)){
      vector<double> log_denom_vec;
    
      for(size_t j = 0; j < distros.size(); j++){
	log_denom_vec.push_back(log(mixing[j]) +
				distros[j].trunc_log_pdf(i));
	if(!finite(log_denom_vec[j])){
	  log_denom_vec[j] = -std::numeric_limits<double>::max()/
	    static_cast<double>(distros.size()); // arbitrary
	}
      }
    
      const double log_denom = log_sum_log_vec(log_denom_vec, 
                                             log_denom_vec.size());
      assert(finite(log_denom));
  
      for(size_t j = 0; j < distros.size(); j++)
	probs[j][i] = exp(log_denom_vec[j] - log_denom);
    }
  }
}


void
ZTNBD_mixture::trunc_maximization_step(const vector<size_t> &vals_hist,
				       const vector< vector<double> > &probs){

  for(size_t i = 0; i < distros.size(); i++){
    distros[i].estim_params(vals_hist, probs[i]);
  }
}

void
ZTNBD_mixture::compute_mixing_w_zeros(const size_t values_size,
				      vector<double> &mixing_w_zeros){
  vector<double> mus;
  vector<double> alphas;
  for(size_t i = 0; i < distros.size(); i++){
    mus.push_back(distros[i].get_mu());
    alphas.push_back(distros[i].get_alpha());
  }

  vector<double> log_genome_lengths_vec;
  for(size_t i =0; i < mixing.size(); i++){
    log_genome_lengths_vec.push_back(log(values_size) + log(mixing[i])
				     -log(1-exp(-log(1+alphas[i]*mus[i])/alphas[i])));
  }
  const double log_total_genome_length = 
    log_sum_log_vec(log_genome_lengths_vec,
		    log_genome_lengths_vec.size());

  for(size_t  i = 0; i < mixing_w_zeros.size(); i++)
    mixing_w_zeros[i] = exp(log_genome_lengths_vec[i]
                            -log_total_genome_length);
}

/*
void
ZTNBD_mixture::compute_probs_w_zeros(const vector<double> &mixing_w_zeros,
				     const vector<size_t> &vals_hist,
				     vector< vector<double> > &probs_w_zeros){
  for(size_t i = 0; i < vals_hist.size(); i++){
    if( (i == 0) || (vals_hist[i] > 0) ){
      vector<double> log_denom_vec;
    
      for(size_t j = 0; j < distros.size(); j++){
	log_denom_vec.push_back(log(mixing_w_zeros[j]) +
				distros[j].log_pdf(i));
	if(!finite(log_denom_vec[j])){
	  log_denom_vec[j] = -std::numeric_limits<double>::max()/
	    static_cast<double>(distros.size()); // arbitrary
	}
      }
    
      const double log_denom = log_sum_log_vec(log_denom_vec, 
                                             log_denom_vec.size());
      assert(finite(log_denom));
  
      for(size_t j = 0; j < distros.size(); j++)
	probs_w_zeros[j][i] = exp(log_denom_vec[j] - log_denom);
    }
  }
}
*/

double 
ZTNBD_mixture::EM_resolve_mix_add_zeros(const double &tol,
					const size_t max_iter,
					const bool VERBOSE,
					vector<size_t> &vals_hist){
  const size_t number_states = distros.size();
  double probs_starting_val = 1/static_cast<double>(number_states);
  vector< vector<double> > probs(number_states, 
                                 vector<double>(vals_hist.size(),
                                                probs_starting_val));

  double score = 0.0;
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iter; ++i){
    trunc_expectation_step(vals_hist, probs);
 
    double expected_zeros = 0.0;
    for(size_t j = 0; j < distros.size(); j++){
      double psuedo_size = 0.0;
      for(size_t l = 1; l < vals_hist.size(); l++){
	psuedo_size += vals_hist[l]*probs[j][l];
      }
      expected_zeros += distros[j].expected_zeros(psuedo_size);
    }
    vals_hist[0] = round(expected_zeros);
    trunc_maximization_step(vals_hist, probs);

    vals_hist[0] = 0;

    trunc_calculate_mixing(vals_hist, probs);

    if(VERBOSE){
      for(size_t j = 0; j < mixing.size(); j++){
	cerr << "params  = " << mixing[j] << "\t" << distros[j].get_mu() 
	     << "\t" << distros[j].get_alpha() << "\n";
      }
      cerr << "\n";
      cerr << "zero probs = ";
      for(size_t j = 0; j < distros.size(); j++){
        cerr << probs[j][0] << ", ";
      }
      cerr << "\n";
    }

    score = trunc_log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;
  }
  const size_t vals_size = accumulate(vals_hist.begin(), 
				      vals_hist.end(), 0);
  double expected_zeros = 0.0;
  for(size_t j = 0; j < distros.size(); j++){
    expected_zeros += distros[j].expected_zeros(mixing[j]*vals_size);
  }
  /*
  vector<double> mixing_w_zeros(mixing);
  compute_mixing_w_zeros(vals_size, mixing_w_zeros);
  */
  compute_Fisher_info(probs, vals_hist, 
		      mixing, expected_zeros);

  return(trunc_log_L(vals_hist));
}

double
ZTNBD_mixture::expected_inverse_sum(const double mean,
                                    const size_t sample_size,
                                    const size_t sum){
  vector<ZTNBD> distro = get_distros();
  vector<double> alphas;
  vector<double> mus;
  vector<double> mixings = get_mixing();
  for(size_t i = 0; i < distro.size(); i++){
    alphas.push_back(distro[i].get_alpha());
    mus.push_back(distro[i].get_mu());
  }
  double expected_MN = 0.0;
  for(size_t i = 0; i< mixings.size(); i++){
    expected_MN += 
      exp(log(mixings[i]) 
	  + log(distro[i].expected_inverse_sum(mean, sample_size,sum)));
  }
  return(expected_MN);
}
double
ZTNBD_mixture::expected_population_size(const size_t sample_size){
  vector<ZTNBD> distro = get_distros();
  vector<double> alphas;
  vector<double> mus;
  vector<double> thetas = get_mixing();
  for(size_t i = 0; i < distro.size(); i++){
    alphas.push_back(distro[i].get_alpha());
    mus.push_back(distro[i].get_mu());
  }
  double expected_pop = 0.0;
  for(size_t i = 0; i < distro.size(); i++){
    expected_pop += 
    thetas[i]*sample_size/(1-exp(-log(1+alphas[i]*mus[i])/alphas[i]));
  }
  return(expected_pop);
}


