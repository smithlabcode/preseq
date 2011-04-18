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

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
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
  const double k = 1/alpha;
  const double w = k/(mu+k);
  double holding_val = 0.0;
  if(val > 0){
    for(size_t j =0; j < val; ++j){
      holding_val += log(k+j);
    }
    holding_val += gsl_sf_lngamma(val);
  }
  return(holding_val + k*log(w) + val*log(1-w));
}

double
NBD::log_L(const vector<size_t> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));
  double log_L = 0;
  const double k = 1/alpha;
  const double w = k/(mu+k);
  for(size_t i = 1; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += log(k+j);

    log_L += vals_hist[i]*gsl_sf_lngamma(i) 
                        + vals_hist[i]*holding_val;

  }
 
  log_L += vals_size*k*log(w);
    
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
  const double k = 1/get_alpha();
  const double w = k/(get_mu()+k);
  double holding_val = 0.0;
  if(val > 0){
    for(size_t j =0; j < val; ++j){
      holding_val += log(k+j);
    }
    holding_val += gsl_sf_lngamma(val);
  }
  return(holding_val + k*log(w) - log(1-pow(w,k)) + val*log(1-w));
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
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));
  double log_L = 0;
  const double k = 1/get_alpha();
  const double w = k/(get_mu()+k);
  for(size_t i = 1; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += log(k+j);

    log_L += vals_hist[i]*gsl_sf_lngamma(i) 
                        + vals_hist[i]*holding_val;

  }
 
  log_L += vals_size*k*log(w) - vals_size*log(1-pow(w,k));
    
  return(log_L);
}

double 
ZTNBD::trunc_pval(const size_t val, const double tol){
  double error = numeric_limits<double>::max();
  double pval = 0.0;
  size_t X = val;
  while(error >= tol){
    error = exp(trunc_log_pdf(X));
    pval += error;
    X++;
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
ZTNBD::expected_inverse_sum(const size_t sample_size,
			    const size_t sum){
  const double alph = get_alpha();
  const double m = get_mu();
  const double prob_zero = pow(1+alph*m,-1/alph);
  const double mean = m/(1-prob_zero);
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
 
    double expected_zeros = 0.0;
    for(size_t j = 0; j < distros.size(); j++){
      double psuedo_size = 0.0;
      for(size_t l = 1; l < vals_hist.size(); l++){
	psuedo_size += vals_hist[l]*probs[j][l];
      }
     }

    maximization_step(vals_hist, probs);
    calculate_mixing(vals_hist, probs);

    score = log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;

  }

  return(log_L(vals_hist));
}




double
obs_fish_info_2nd_theta(const vector< vector<double> > &probs,
			const vector<size_t> &vals_hist,
			const vector<double> &mixing,
			const size_t indx,
			const double expected_zeros){
  vector<double> log_pos_vec;
  vector<double> log_neg_vec;
  /*
  if(finite(log(probs[indx][0])) && expected_zeros > 0){
    log_pos_vec.push_back(log(expected_zeros) + 2*log(probs[indx][0])
			    -2*log(mixing[indx]));
  }
  if(finite(log(probs[probs.size()-1][0])) && expected_zeros > 0){
    log_pos_vec.push_back(2*log(probs[probs.size()-1][0])
			  +log(expected_zeros)
			  -2*log(mixing.back()));
  }
  if(finite(log(probs[indx][0])) &&
     finite(log(probs[probs.size()-1][0])) && expected_zeros > 0){
    log_neg_vec.push_back(log(probs[indx][0])+log(probs[probs.size()-1][0])
			  +log(expected_zeros)-log(mixing[indx])
			  -log(mixing.back()));
  }
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      if(finite(log(probs[indx][i]))){
	log_pos_vec.push_back(log(vals_hist[i]) + 2*log(probs[indx][i])
				-2*log(mixing[indx]));
      }
      if(finite(log(probs[probs.size()-1][i]))){
	log_pos_vec.push_back(log(vals_hist[i])
				+2*log(probs[probs.size()-1][i])
				-2*log(mixing.back()));
      }
      if(finite(log(probs[indx][i])) && finite(log(probs[probs.size()-1][i]))){
	log_neg_vec.push_back(log(vals_hist[i])+log(probs[indx][i])
			      +log(probs[probs.size()-1][i])
			      -log(mixing[indx])-log(mixing.back()));
      }
    }
  }
  return(exp(log_sum_log_vec(log_pos_vec, log_pos_vec.size()))
	 -exp(log_sum_log_vec(log_neg_vec, log_neg_vec.size())));
}

double 
obs_fish_info_mixed_theta(const vector< vector<double> > &probs,
			  const vector<size_t> &vals_hist,
			  const vector<double> &mixing,
			  const double expected_zeros,
			  const size_t indx1,
			  const size_t indx2){
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  /*
  if(expected_zeros > 0){
    if(finite(log(probs[indx1][0])) && finite(log(probs[indx2][0]))){
      pos_log_vec.push_back(log(expected_zeros) + log(probs[indx1][0])
			    +log(probs[indx2][0]) -log(mixing[indx1])
			    -log(mixing[indx2]));
    }
    if(finite(log(probs[probs.size()-1][0]))){
      pos_log_vec.push_back(log(expected_zeros)
			    +2*log(probs[probs.size()-1][0])
			    -2*log(mixing.back()));
    }
    if(finite(log(probs[indx1][0])) && finite(log(probs[indx2][0]))
       && finite(log(probs[probs.size()-1][0]))){
      neg_log_vec.push_back(log(expected_zeros)
			    +log(probs[probs.size()-1][0])
			    -log(mixing.back())
			    +log(probs[indx1][0]/mixing[indx1]
				 +probs[indx2][0]/mixing[indx2]));
    }
  }
  */
  for(size_t i =1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_prob_1 = log(probs[indx1][i]);
      const double log_prob_2 = log(probs[indx2][i]);
      const double log_last_prob = log(probs[probs.size()-1][i]);
      if(finite(log_prob_1) &&
	 finite(log_prob_2) &&
	 finite(log_last_prob)){
	neg_log_vec.push_back(log_vals_hist
			      +log_last_prob
			      -log(mixing.back())
			      +log(probs[indx1][i]/mixing[indx1]
				   +probs[indx2][i]/mixing[indx2]));
      }
      if(finite(log_prob_1) &&
	 finite(log_prob_2)){
	pos_log_vec.push_back(log_vals_hist
			      +log_prob_1
			      +log_prob_2
			      -log(mixing[indx1])-log(mixing[indx2]));
      }
      if(finite(log_last_prob)){
	pos_log_vec.push_back(log_vals_hist
			      +2*log_last_prob
			      -2*log(mixing.back()));
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double
obs_fish_info_2nd_mu(const vector< vector<double> > &probs,
		     const vector<size_t> &vals_hist,
		     const vector<ZTNBD> &distros,
		     const double expected_zeros,
		     const size_t indx){
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();
  const double plus_alpha_mu = 1+alpha*mu;
  /*
  if(expected_zeros > 0 && finite(log(probs[indx][0]))){
    neg_log_vec.push_back(log(expected_zeros)+log(probs[indx][0])
			  +log(alpha)-2*log(plus_alpha_mu));
    neg_log_vec.push_back(log(expected_zeros)+log(probs[indx][0])
			  -2*log(plus_alpha_mu));
    pos_log_vec.push_back(log(expected_zeros)+2*log(probs[indx][0])
			  -2*log(plus_alpha_mu));
  }
  */
  for(size_t i =1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_prob = log(probs[indx][i]);
      if(finite(log_prob)){
	pos_log_vec.push_back(log_vals_hist+log_prob+log(i)
			      -2*log(mu));
	neg_log_vec.push_back(log_vals_hist+log_prob+log(alpha)
			      +log(1+alpha*i)-2*log(plus_alpha_mu));
	neg_log_vec.push_back(log_vals_hist+log_prob+2*log(i)
			      -2*log(mu));
	pos_log_vec.push_back(log_vals_hist+log_prob+log(2)
			      +log(i)+log(1+alpha*i)
			      -log(mu)-log(plus_alpha_mu));
	neg_log_vec.push_back(log_vals_hist+log_prob
			      +2*log(1+alpha*i)-2*log(plus_alpha_mu));
	pos_log_vec.push_back(log_vals_hist+2*log_prob
			      +2*log(i)-2*log(mu));
	pos_log_vec.push_back(log_vals_hist+2*log_prob
			      +2*log(1+alpha*i)
			      -2*log(plus_alpha_mu));
	neg_log_vec.push_back(log_vals_hist+2*log_prob+log(2)
			      +log(i)+log(1+alpha*i)-log(mu)
			      -log(plus_alpha_mu));
      
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double 
obs_fish_info_mixed_mu(const vector< vector<double> > &probs,
		       const vector<size_t> &vals_hist,
		       const vector<ZTNBD> &distros,
		       const double expected_zeros,
		       const size_t indx1,
		       const size_t indx2){
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  const double mu1 = distros[indx1].get_mu();
  const double mu2 = distros[indx2].get_mu();
  const double alpha1 = distros[indx1].get_alpha();
  const double alpha2 = distros[indx2].get_alpha();
  /*
  if(finite(log(probs[indx1][0])) && finite(log(probs[indx2][0]))){
    pos_log_vec.push_back(log(expected_zeros)+log(probs[indx1][0])
			  +log(probs[indx2][0])-log(1+alpha1*mu1)
			  -log(1+alpha1*mu1));
  }
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_probs1 = log(probs[indx1][i]);
      const double log_probs2 = log(probs[indx2][i]);
      const double log_vals_hist = log(vals_hist[i]);
      if(finite(log_probs1) && finite(log_probs2)){
	pos_log_vec.push_back(log_vals_hist+log_probs1+log_probs2
			      +2*log(i)-log(mu1)-log(mu2));
	pos_log_vec.push_back(log_vals_hist+log_probs1+log_probs2
			      +log(1+alpha1*i)+log(1+alpha2*i)
			      -log(1+alpha1*mu1)-log(1+alpha2*mu2));
	neg_log_vec.push_back(log_vals_hist+log_probs1+log_probs2
			      +log(i)-log(mu1)+log(1+alpha2*i)
			      -log(1+alpha2*mu2));
	neg_log_vec.push_back(log_vals_hist+log_probs1+log_probs2
			      +log(i)-log(mu2)+log(1+alpha1*i)
			      -log(1+alpha1*mu1));
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));

}

double
obs_fish_info_2nd_alpha(const vector< vector<double> > &probs,
			const vector<size_t> &vals_hist,
			const vector<ZTNBD> &distros,
			const double expected_zeros,
			const size_t indx){
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;

  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();
  for(size_t i = 1; i < vals_hist.size(); i++){
    const double log_prob = log(probs[indx][i]);
    const double log_1_minus_prob = log(1-probs[indx][i]);

    /*
    if(i ==0 && finite(log_prob)){
      if(log(1+alpha*mu) > 0){
	pos_log_vec.push_back(log(expected_zeros) + log_prob + 
			      log(2) + log(log(1+alpha*mu)) - 3*log(alpha));
	if(finite(log_1_minus_prob)){
	  neg_log_vec.push_back(log(expected_zeros) + log_prob
				+log_1_minus_prob
				+2*log(log(1+alpha*mu))-4*log(alpha));
	  pos_log_vec.push_back(log(expected_zeros)+log_prob
				+log_1_minus_prob
				+log(log(1+alpha*mu))+log(mu)
				-3*log(alpha)-log(1+alpha*mu));
	}
      }
      else if(log(1+alpha*mu) < 0){
	neg_log_vec.push_back(log(expected_zeros)+log_prob+log(2)
			      + log(log(1+alpha*mu)) - 3*log(alpha));
	if(finite(log_1_minus_prob)){
	  pos_log_vec.push_back(log(expected_zeros) + log_prob
				+log_1_minus_prob
				+2*log(log(1+alpha*mu))-4*log(alpha)
				+2*log(log(1+alpha*mu))-4*log(alpha));
	  neg_log_vec.push_back(log(expected_zeros)+log_prob
				+log_1_minus_prob
				+log(log(1+alpha*mu))+log(mu)
				-3*log(alpha)-log(1+alpha*mu));
	}
      }
      neg_log_vec.push_back(log(expected_zeros)+log_prob+log(2)
			    +log(mu)-2*log(alpha)-log(1+alpha*mu)); 
      if(finite(log_1_minus_prob)){
	neg_log_vec.push_back(log(expected_zeros)+log_prob+log_1_minus_prob
			      +2*log(mu)-2*log(alpha)-2*log(1+alpha*mu));
      }
    }
    
  //not zero case
    else{
    */
      if(vals_hist[i] > 0){

	const double log_vals_hist = log(vals_hist[i]);
	double first_deriv_sum = 0.0;
	double second_deriv_sum = 0.0;
	for(size_t j = 0; j < i; j++){
	  double holding_val = j/(1+alpha*j);
	  first_deriv_sum += holding_val;
	  second_deriv_sum += holding_val*holding_val;
	}

	if(finite(log_prob)){
	  pos_log_vec.push_back(log_vals_hist+log_prob+
				log(second_deriv_sum));
	  pos_log_vec.push_back(log_vals_hist+log_prob
				+log(2)-3*log(alpha)+
				log(log(1+alpha*mu)));
	  neg_log_vec.push_back(log_vals_hist+log_prob+log(2)
				+log(mu)-2*log(alpha)-log(1+alpha*mu));
	  neg_log_vec.push_back(log_vals_hist+log_prob+2*log(mu)
				+log(1+alpha*i)-log(alpha)
				-2*log(1+alpha*mu));
	  pos_log_vec.push_back(log_vals_hist+log_prob+log(i)
				-2*log(alpha));
	  if(finite(log_1_minus_prob)){
	    double first_deriv = first_deriv_sum + exp(log(i)-log(alpha))
	      +exp(log(log(1+alpha*mu))-2*log(alpha))-
	      exp(log(mu)+log(1+alpha*i)-log(alpha)-log(1+alpha*mu));
	    double first_deriv_2 = first_deriv*first_deriv;
	    if(first_deriv_2 > 0){
	      neg_log_vec.push_back(log_vals_hist+log_prob
				    +log_1_minus_prob + log(first_deriv_2)); 
	    }
	  }
	}
      }
      // }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double
obs_fish_info_mixed_alpha(const vector< vector<double> > &probs,
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
  /*
  terms_vec.push_back(exp(log(expected_zeros)+log(probs[indx1][0])
			 +log(probs[indx2][0]))*
		     (log(1+alpha1*mu1)/(alpha1*alpha1)
		      -mu1/(alpha1*(1+alpha1*mu1)))*
		     (log(1+alpha2*mu2)/(alpha2*alpha2)
		      -mu2/(alpha2*(1+alpha2*mu2))));
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      double first_deriv_sum1 = 0.0;
      double first_deriv_sum2 = 0.0;
      for(size_t j = 0; j < i; j++){
	first_deriv_sum1 += j/(1+alpha1*j);
	first_deriv_sum2 += j/(1+alpha2*j);
      }
      terms_vec.push_back(exp(log(vals_hist[i])+log(probs[indx1][i])
			      +log(probs[indx2][i]))*
			  (first_deriv_sum1 + log(1+alpha1*mu1)/(alpha1*alpha1)
			   +i/alpha1-mu1*(1+alpha1*i)/(alpha1*(1+alpha1*mu1)))*
			  (first_deriv_sum2+log(1+alpha2*mu2)/(alpha2*alpha2)
			   +i/alpha2-mu2*(1+alpha2*i)/(alpha2*(1+alpha2*mu2))));
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}


double
obs_fish_info_dtheta_dmu_same_indx(const vector< vector<double> > &probs,
				   const vector<size_t> &vals_hist,
				   const vector<ZTNBD> &distros,
				   const vector<double> &mixing,
				   const double expected_zeros,
				   const size_t indx){
  const double last_mix = mixing.back();
  const double mix = mixing[indx];
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  /*
  if(finite(log(probs[indx][0])) && expected_zeros > 0){
    pos_log_vec.push_back(log(expected_zeros)+log(probs[indx][0])
			  -log(1+alpha*mu));
    neg_log_vec.push_back(log(expected_zeros)+2*log(probs[indx][0])
			  -log(1+alpha*mu));
    if(finite(log(probs[probs.size()-1][0]))){
      pos_log_vec.push_back(log(expected_zeros)+log(probs[indx][0])
			    +log(probs[probs.size()-1][0])
			    -log(last_mix)-log(1+alpha*mu));
    }
  }
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_prob = log(probs[indx][i]);
      const double log_last_prob = log(probs[probs.size()-1][i]);
      const double log_vals_hist = log(vals_hist[i]);
      if(finite(log_prob)){
	neg_log_vec.push_back(log_vals_hist+log_prob-log(mix)
			      +log(i)-log(mu));
	pos_log_vec.push_back(log_vals_hist+log_prob+log(i)
			      +log(alpha)-log(mix)
			      -log(1+alpha*mu));
	pos_log_vec.push_back(log_vals_hist+log_prob-log(mix)
			      -log(1+alpha*mu));
	pos_log_vec.push_back(log_vals_hist+2*log_prob+log(i)
			      -log(mix)-log(mu));
	neg_log_vec.push_back(log_vals_hist+2*log_prob+log(i)
			      +log(alpha)-log(mix)
			      -log(1+alpha*mu));
	neg_log_vec.push_back(log_vals_hist+2*log_prob-log(mix)
			      -log(1+alpha*mu));
	if(finite(log_last_prob)){
	  neg_log_vec.push_back(log_vals_hist+log_prob+log_last_prob
				-log(last_mix)+log(i)-log(mu));
	  pos_log_vec.push_back(log_vals_hist+log_prob+log_last_prob
				-log(last_mix)+log(i*alpha+1)
				-log(1+alpha*mu));
	}
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double
obs_fish_info_dtheta_dmu_diff_indx(const vector< vector<double> > &probs,
				   const vector<size_t> &vals_hist,
				   const vector<ZTNBD> &distros,
				   const vector<double> &mixing,
				   const double expected_zeros,
				   const size_t indx_mix,
				   const size_t indx_mu){
  const double mix = mixing[indx_mix];
  const double last_mix = mixing.back();
  const double mu = distros[indx_mu].get_mu();
  const double alpha = distros[indx_mu].get_alpha();

  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  /*
  if(finite(log(probs[indx_mix][0])) && finite(log(probs[indx_mu][0]))
     && expected_zeros > 0){
    neg_log_vec.push_back(log(expected_zeros)+log(probs[indx_mix][0])
			  +log(probs[indx_mu][0])-log(mix)
			  -log(1+alpha*mu));
  }
  if(finite(log(probs[probs.size()-1][0])) && finite(log(probs[indx_mu][0]))
     && expected_zeros > 0){
    pos_log_vec.push_back(log(expected_zeros) + log(probs[probs.size()-1][0])
			  +log(probs[indx_mu][0])-log(last_mix)
			  -log(1+alpha*mu));
  }
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_prob_mix = log(probs[indx_mix][i]);
      const double log_prob_mu = log(probs[indx_mu][i]);
      const double log_prob_last = log(probs[probs.size()-1][i]);
      const double log_vals_hist = log(vals_hist[i]);
      if(finite(log_prob_mix) && finite(log_prob_mu)){
	pos_log_vec.push_back(log_vals_hist+log_prob_mix+log_prob_mu
			      -log(mix)+log(i) -log(mu));
	neg_log_vec.push_back(log_vals_hist+log_prob_mix+log_prob_mu
			      -log(mix)+log(1+i*alpha)
			      -log(1+alpha*mu));
      }
      if(finite(log_prob_last) && finite(log_prob_mu)){
	pos_log_vec.push_back(log_vals_hist+log_prob_last+log_prob_mu
			      -log(last_mix)+log(i) -log(mu));
	neg_log_vec.push_back(log_vals_hist+log_prob_last+log_prob_mu
			      -log(last_mix)+log(1+i*alpha)
			      -log(1+alpha*mu));
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double
obs_fish_info_dtheta_d_last_mu(const vector< vector<double> > &probs,
			       const vector<size_t> &vals_hist,
			       const vector<ZTNBD> &distros,
			       const vector<double> &mixing,
			       const double expected_zeros,
			       const size_t indx){
  const double mu = distros[distros.size()-1].get_mu();
  const double alpha = distros[distros.size()-1].get_alpha();
  const double last_mix = mixing.back();
  const double mix = mixing[indx];

  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  /*
  if(finite(log(probs[probs.size()-1][0]))){
    neg_log_vec.push_back(log(expected_zeros)
			  +log(probs[probs.size()-1][0])
			  -log(last_mix)-log(1+alpha*mu));
    pos_log_vec.push_back(log(expected_zeros)
			  +2*log(probs[probs.size()-1][0])
			  -log(last_mix)-log(1+alpha*mu));
    if(finite(log(probs[indx][0]))){
      neg_log_vec.push_back(log(expected_zeros)
			    +log(probs[indx][0])
			    +log(probs[probs.size()-1][0])
			    -log(mix)-log(1+alpha*mu));
    }
  }
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_last_prob = log(probs[probs.size()-1][i]);
      const double log_prob = log(probs[indx][i]);
      const double log_vals_hist = log(vals_hist[i]);
      if(finite(log_last_prob)){
	pos_log_vec.push_back(log_vals_hist+log_last_prob
			      -log(last_mix)+log(i)-log(mu));
	neg_log_vec.push_back(log_vals_hist+log_last_prob
			      -log(last_mix)+log(1+alpha*i)
			      -log(1+alpha*mu));
	neg_log_vec.push_back(log_vals_hist+2*log_last_prob
			      -log(last_mix)+log(i)-log(mu));
	pos_log_vec.push_back(log_vals_hist+2*log_last_prob
			      -log(last_mix)+log(1+alpha*i)
			      -log(1+alpha*mu));
	if(finite(log_prob)){
	  pos_log_vec.push_back(log_vals_hist+log_prob+log_last_prob
				-log(mix)+log(i)-log(mu));
	  neg_log_vec.push_back(log_vals_hist+log_prob+log_last_prob
				-log(mix)+log(1+alpha*i)
				-log(1+alpha*mu));
	}
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}
      
double
obs_fish_info_dtheta_dalpha_same_indx(const vector< vector<double> > &probs,
				      const vector<size_t> &vals_hist,
				      const vector<ZTNBD> &distros,
				      const vector<double> &mixing,
				      const double expected_zeros,
				      const size_t indx){
  const double mu = distros[indx].get_mu();
  const double alpha = distros[indx].get_alpha();
  const double mix = mixing[indx];
  const double last_mix = mixing.back();
  vector<double> terms_vec;
  /*
  terms_vec.push_back(-expected_zeros*probs[indx][0]*(1-probs[indx][0])*
		      (log(1+alpha*mu)/(alpha*alpha)-mu/(alpha*(1+alpha*mu)))
		      /mix);
  terms_vec.push_back(-expected_zeros*probs[indx][0]*probs[probs.size()-1][0]*
		      log(1+alpha*mu)/(alpha*alpha)-
		      mu/(alpha*(1+alpha*mu))/last_mix);
  */
  for(size_t i  = 0; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      double inner_sum = 0;
      for(size_t j = 0; j < i; j++){
	inner_sum += j/(1+alpha*j);
      }
      const double prob = probs[indx][i];
      const double last_prob = probs[probs.size()-1][i];
      terms_vec.push_back(-vals_hist[i]*prob*(1-prob)*
			  (inner_sum+log(1+alpha*mu)/(alpha*alpha)+i/alpha
			   -mu*(1+alpha*i)/(alpha*(1+alpha*mu)))/mix);
      terms_vec.push_back(-vals_hist[i]*prob*last_prob*
			  (inner_sum+log(1+alpha*mu)/(alpha*alpha)+i/alpha
			   -mu*(1+alpha*i)/(alpha*(1+alpha*mu)))/last_mix);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
obs_fish_info_dtheta_dalpha_diff_indx(const vector< vector<double> > &probs,
				      const vector<size_t> &vals_hist,
				      const vector<ZTNBD> &distros,
				      const vector<double> &mixing,
				      const double expected_zeros,
				      const size_t indx_mix,
				      const size_t indx_alpha){
  const double alpha = distros[indx_alpha].get_alpha();
  const double mu = distros[indx_alpha].get_mu();
  const double mix = mixing[indx_mix];
  const double last_mix = mixing.back();

  vector<double> terms_vec;
  /*
  terms_vec.push_back(expected_zeros*probs[indx_mix][0]*probs[indx_alpha][0]*
		      (log(1+alpha*mu)/(alpha*alpha)
		       -mu/(alpha*(1+alpha*mu)))/mix);
  terms_vec.push_back(-expected_zeros*probs[indx_alpha][0]*
		      probs[probs.size()-1][0]*
		      (log(1+alpha*mu)/(alpha*alpha)
		       -mu/(alpha*(1+alpha*mu)))/last_mix);
  */
  for(size_t  i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double prob_alpha = probs[indx_alpha][i];
      const double prob_mix = probs[indx_mix][i];
      const double prob_last = probs[probs.size()-1][i];
      double inner_sum =  0.0;
      for(size_t j =0; j < i; j++){
	inner_sum += j/(1+alpha*j);
      }
      const double hist = static_cast<double>(vals_hist[i]);
      terms_vec.push_back(hist*prob_mix*prob_alpha*
			  (inner_sum+ log(1+alpha*mu)/(alpha*alpha)+i/alpha
			   -mu*(1+alpha*i)/(alpha*(1+alpha*mu)))/mix);
      terms_vec.push_back(-hist*prob_last*prob_alpha*
			  (inner_sum+ log(1+alpha*mu)/(alpha*alpha)+i/alpha
			   -mu*(1+alpha*i)/(alpha*(1+alpha*mu)))/last_mix);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double 
obs_fish_info_dtheta_d_last_alpha(const vector< vector<double> > &probs,
				  const vector<size_t> &vals_hist,
				  const vector<ZTNBD> &distros,
				  const vector<double> &mixing,
				  const double expected_zeros,
				  const size_t indx){
  const double alpha = distros[distros.size()-1].get_alpha();
  const double mu = distros[distros.size()-1].get_mu();
  const double last_mix = mixing.back();
  const double mix = mixing[indx];
  double last_prob = probs[probs.size()-1][0];
  vector<double> terms_vec;
  /*
  terms_vec.push_back(expected_zeros*last_prob*(1-last_prob)*
		      (log(1+alpha*mu)/(alpha*alpha)
		       -mu/(alpha*(1+alpha*mu)))/last_mix);
  terms_vec.push_back(expected_zeros*last_prob*probs[indx][0]*
		      (log(1+alpha*mu)/(alpha*alpha)
		       -mu/(alpha*(1+alpha*mu)))/mix);
  */
  for(size_t  i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      last_prob = probs[probs.size()-1][i];
      const double prob = probs[indx][i];
      const double hist = static_cast<double>(vals_hist[i]);
      double inner_sum = 0.0;
      for(size_t j = 0; j < i; j++){
	inner_sum += j/(1+alpha*j);
      }
      terms_vec.push_back(hist*last_prob*(1-last_prob)*
			 (inner_sum+log(1+alpha*mu)/(alpha*alpha)+i/alpha
			  -mu*(1+alpha*i)/(alpha*(1+alpha*mu)))/last_mix);
      terms_vec.push_back(hist*last_prob*prob*
			 (inner_sum+log(1+alpha*mu)/(alpha*alpha)+i/alpha
			  -mu*(1+alpha*i)/(alpha*(1+alpha*mu)))/mix);
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
obs_fish_info_dmu_dalpha_same_indx(const vector< vector<double> > &probs,
				   const vector<size_t> &vals_hist,
				   const vector<ZTNBD> &distros,
				   const double expected_zeros,
				   const size_t indx){
  const double alpha = distros[indx].get_alpha();
  const double mu = distros[indx].get_mu();
  const double log_mu=log(mu);
  const double log_alpha = log(alpha);
  const double log_alpha_mu = log(1+alpha*mu);

  vector<double> terms_vec;

  const double log_prob_zero = log(probs[indx][0]);
  const double log_1_minus_prob_zero = log(1-probs[indx][0]);
  const double log_e_zeros = log(expected_zeros);
  /*
  if(finite(log_prob_zero) && expected_zeros > 0){
    terms_vec.push_back(-exp(log_e_zeros + log_prob_zero + log_mu
			     -2*log_alpha_mu));
    terms_vec.push_back(exp(log_e_zeros+log_prob_zero
			    +log_1_minus_prob_zero-log_alpha_mu
			    -2*log_alpha));
    terms_vec.push_back(-exp(log_e_zeros+log_prob_zero
			     +log_1_minus_prob_zero-log_alpha_mu
			     +log_mu-log_alpha-log_alpha_mu));
  }
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_prob = log(probs[indx][i]);
      const double log_one_minus_prob = log(1-probs[indx][i]);
      const double log_alpha_i = log(alpha*i+1);
      const double log_vals_hist = log(vals_hist[i]);
      if(finite(log_prob)){
	terms_vec.push_back(exp(log_vals_hist+log_prob+log(i)
				-log_alpha_mu));
	terms_vec.push_back(-exp(log_vals_hist+log_prob+log(mu)
				 +log_alpha_i-2*log_alpha_mu));
	if(finite(log_one_minus_prob)){
	  double inner_sum = 0.0;
	  for(size_t j = 0; j < i; j++)
	    inner_sum += j/(1+j*alpha);
	  const double alpha_deriv = inner_sum + i/alpha
	    +log(1+alpha*mu)/(alpha*alpha)
	    -mu*(1+alpha*i)/(alpha*(1+alpha*mu));
	  const double mu_deriv = i/mu - (1+i*alpha)/(1+alpha*mu);
	  terms_vec.push_back(-exp(log_vals_hist+log_prob
				   +log_one_minus_prob)*alpha_deriv*mu_deriv);
	}
      }
    }
  }
  return(accumulate(terms_vec.begin(), terms_vec.end(), 0.0));
}

double
obs_fish_info_dmu_dalpha_diff_indx(const vector< vector<double> > &probs,
				   const vector<size_t> &vals_hist,
				   const vector<ZTNBD> &distros,
				   const double expected_zeros,
				   const size_t indx_mu,
				   const size_t indx_alpha){

  const double alpha1 = distros[indx_mu].get_alpha();
  const double mu1 = distros[indx_mu].get_mu();
  const double alpha2 = distros[indx_alpha].get_alpha();
  const double mu2 = distros[indx_alpha].get_mu();
  const double log_alpha_mu_1 = log(1+alpha1*mu1);
  const double log_alpha_mu_2 = log(1+alpha2*mu2);

  vector<double> terms_vec;
  /*
  terms_vec.push_back(exp(log(expected_zeros)+log(probs[indx_mu][0])
			  +log(probs[indx_alpha][0]) - log_alpha_mu_1)*
		      (-log_alpha_mu_2/(alpha2*alpha2)+
		       exp(log(mu2)-log(alpha2)-log_alpha_mu_2)));
  */
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_probs = 
	log(probs[indx_mu][i])+log(probs[indx_alpha][i]);
      double inner_sum = 0.0;
      for(size_t j = 0; j < i; j++)
	inner_sum += j/(1+alpha2*j);
      
      if(finite(log_probs)){
	const double alpha_deriv = inner_sum + i/alpha2
	    +log(1+alpha2*mu2)/(alpha2*alpha2)
	    -mu2*(1+alpha2*i)/(alpha2*(1+alpha2*mu2));
	  const double mu_deriv = i/mu1 - (1+i*alpha1)/(1+alpha1*mu1);
	  terms_vec.push_back(-exp(log_vals_hist+log_probs)
			      *alpha_deriv*mu_deriv);
      }
    }
  }
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
    for(size_t j = i; j < num_mixs-1; j++){
      if(j == i){
	holding_fish[i][j] = obs_fish_info_2nd_theta(probs, vals_hist,
						     mixing, i, 
						     expected_zeros);
	
      }
      else{
	holding_fish[i][j] = obs_fish_info_mixed_theta(probs, vals_hist,
						       mixing, expected_zeros,
						       i, j);
	holding_fish[j][i] = holding_fish[i][j];
      }
    }
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = i; j < num_mixs; j++){
      if(j==i){
	holding_fish[i+num_mixs-1][j+num_mixs-1] = 
	  obs_fish_info_2nd_mu(probs, vals_hist, distros,
			       expected_zeros, i);

      }
      else{
	holding_fish[i+num_mixs-1][j+num_mixs-1] = 
	  obs_fish_info_mixed_mu(probs, vals_hist, distros,
				 expected_zeros, i, j);
	holding_fish[j+num_mixs-1][i+num_mixs-1] = 
	  holding_fish[i+num_mixs-1][j+num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = i; j < num_mixs; j++){
      if(j==i){
	holding_fish[i+2*num_mixs-1][j+2*num_mixs-1] = 
	  obs_fish_info_2nd_alpha(probs, vals_hist, distros,
				  expected_zeros, i);
      }
      else{
	holding_fish[i+2*num_mixs-1][j+2*num_mixs-1] = 
	  obs_fish_info_mixed_alpha(probs, vals_hist, distros,
				    expected_zeros, i, j);
	holding_fish[j+2*num_mixs-1][i+2*num_mixs-1] = 
	  holding_fish[i+2*num_mixs-1][j+2*num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = 0; j < num_mixs-1; j++){
      if(j==i){
	holding_fish[i][j+num_mixs-1] = 
	  obs_fish_info_dtheta_dmu_same_indx(probs, vals_hist,
					     distros, mixing, 
					     expected_zeros, i);
	holding_fish[j+num_mixs-1][i] = 
	  holding_fish[i][j+num_mixs-1];
      }
      else{
	holding_fish[i][j+num_mixs-1] = 
	  obs_fish_info_dtheta_dmu_diff_indx(probs, vals_hist,
					     distros, mixing, 
					     expected_zeros, i, j);
	holding_fish[j+num_mixs-1][i] = 
	  holding_fish[i][j+num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    holding_fish[i][2*num_mixs-2] = 
      obs_fish_info_dtheta_d_last_mu(probs, vals_hist, distros,
				     mixing, expected_zeros, i);
    holding_fish[2*num_mixs-2][i] = 
      holding_fish[i][2*num_mixs-2];
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    for(size_t j = 0; j < num_mixs-1; j++){
      if(j==i){
	holding_fish[i][j+2*num_mixs-1] = 
	  obs_fish_info_dtheta_dalpha_same_indx(probs, vals_hist,
						distros, mixing,
						expected_zeros, i);
	holding_fish[j+2*num_mixs-1][i] = 
	  holding_fish[i][j+2*num_mixs-1];
      }
      else{
	holding_fish[i][j+2*num_mixs-1] = 
	  obs_fish_info_dtheta_dalpha_diff_indx(probs, vals_hist,
						distros, mixing, 
						expected_zeros, i, j);
	holding_fish[j+2*num_mixs-1][i] = 
	  holding_fish[i][j+2*num_mixs-1];
      }
    }
  }
  for(size_t i = 0; i < num_mixs-1; i++){
    holding_fish[i][3*num_mixs-2] = 
      obs_fish_info_dtheta_d_last_alpha(probs, vals_hist,
					distros, mixing, 
					expected_zeros, i);
    holding_fish[3*num_mixs-2][i] = holding_fish[i][3*num_mixs-2];
  }
  for(size_t i = 0; i < num_mixs; i++){
    for(size_t j = i; j < num_mixs; j++){
      if(j==i){
	holding_fish[2*num_mixs-1+i][num_mixs-1+j] = 
	  obs_fish_info_dmu_dalpha_same_indx(probs, vals_hist,
					     distros, expected_zeros, i);
	holding_fish[num_mixs-1+j][2*num_mixs-1+i] = 
	  holding_fish[2*num_mixs-1+i][num_mixs-1+j];
      }
      else{
	holding_fish[2*num_mixs-1+i][num_mixs-1+j] =
	  obs_fish_info_dmu_dalpha_diff_indx(probs, vals_hist,
					     distros, expected_zeros,
					     i, j);
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
ZTNBD_mixture::trunc_expectation_step(const vector<size_t> &vals_hist,
				      vector< vector<double> > &probs){
  for(size_t i = 0; i < vals_hist.size(); i++){
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


void
ZTNBD_mixture::trunc_maximization_step(const vector<size_t> &vals_hist,
				       const vector< vector<double> > &probs){

  for(size_t i = 0; i < distros.size(); i++){
    distros[i].estim_params(vals_hist, probs[i]);
  }
}

void
ZTNBD_mixture::compute_mixing_w_zeros(const size_t sample_size,
				      vector<double> &mixing){
  vector<double> log_genome_lengths;
  vector<ZTNBD> distros= get_distros();
  for(size_t i = 0; i < mixing.size(); i++){
    log_genome_lengths.push_back(log(mixing[i])+log(sample_size)
				 -log(1-exp(-log(1+distros[i].get_mu()*
						 distros[i].get_alpha())
					    /distros[i].get_alpha())));
  }
  double log_total_length = log_sum_log_vec(log_genome_lengths,
					    log_genome_lengths.size());
  for(size_t i = 0; i < mixing.size(); i++){
    mixing[i] = exp(log_genome_lengths[i]-log_total_length);
  }
}

double 
ZTNBD_mixture::EM_resolve_mix_add_zeros(const double &tol,
					const size_t max_iter,
					vector<size_t> &vals_hist){
  const size_t vals_size = accumulate(vals_hist.begin(), vals_hist.end(), 0);
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

    score = trunc_log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;
  }
  //Compute Fisher_info, remember wrt to q, not theta
    vector<double> mixing_w_zeros(mixing);
   /*   compute_mixing_w_zeros(vals_size, mixing_w_zeros); */
  double expected_zeros = 0.0;
  for(size_t j = 0; j < distros.size(); j++){
    double psuedo_size = 0.0;
    for(size_t l = 1; l < vals_hist.size(); l++){
      psuedo_size += vals_hist[l]*probs[j][l];
    }
    expected_zeros += distros[j].expected_zeros(psuedo_size);
  }
  compute_Fisher_info(probs, vals_hist,mixing_w_zeros , expected_zeros);

  return(trunc_log_L(vals_hist));
}

double
ZTNBD_mixture::expected_inverse_sum(const size_t sample_size,
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
	  + log(distro[i].expected_inverse_sum(sample_size,sum)));
  }
  return(expected_MN);
}

