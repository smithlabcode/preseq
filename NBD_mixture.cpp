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
  for(size_t i = 1; i < vals_hist.size(); i++){
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

double 
ZTNBD_mixture::EM_resolve_mix_add_zeros(const double &tol,
					const size_t max_iter,
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

    score = trunc_log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;

  }

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

