/*    ZTNB:
 *
 *    Copyright (C) 2013 University of Southern California and
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
#include "ZTNB.hpp"

#include <smithlab_utils.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>


using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;

using std::max;
using std::fabs;
using std::accumulate;

static double 
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


static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

static inline double
compute_mean(const vector<double> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(),
				   vals_hist.end(), 0));
  double mean = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}


const double ZTNBD::max_allowed_alpha = 1000;
const double ZTNBD::min_allowed_alpha = 1e-20;
const double ZTNBD::tolerance = 1e-10; 

double
ZTNBD::expected_zeros(const double distinct){
  const double alpha = get_alpha();
  const double mu = get_mu();
  const double prob_zero = pow(1+alpha*mu, -1/alpha);
  const double expected_zeros = distinct*(prob_zero/(1-prob_zero));
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


void 
ZTNBD::estim_params(const vector<double> &vals_hist){
  mu = compute_mean(vals_hist);
  //mu= (1/n)sum(x), accumulate takes the sum of vals.begin
  
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), 
				   vals_hist.end(), 0));
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  
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
ZTNBD::estim_params(const vector<double> &vals_hist,
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
  
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(pseudo_hist,mu, a_mid, pseudo_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid; 
  set_helpers();
}


double
ZTNBD::trunc_log_L(const vector<double> &vals_hist){
  double log_L = 0.0;
  const double prob_zero = exp(-log(1+alpha*mu)/alpha);
    double holding_val = 0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    holding_val += log(1 + alpha*(i - 1));
    log_L += 
      vals_hist[i]*(holding_val - gsl_sf_lngamma(i + 1)
		    +i*log(mu) - (i + 1.0/alpha)*log(1 + alpha*mu)
		    - log(1 - prob_zero));
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
		      vector<double> &vals_hist){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(),0));
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();
  double score = 0.0;
  for(size_t i = 0; i < max_iter; i++){
    vals_hist[0] = expected_zeros(vals_size);

    estim_params(vals_hist);

    vals_hist[0] = 0.0;
    score = trunc_log_L(vals_hist);
    error = fabs((score - prev_score)/score);
    if(error < tol)
      break;
    prev_score = score;
  }
  return(trunc_log_L(vals_hist));
}


// for a fixed alpha
double
ZTNBD::EM_estim_mu_fixed_alpha(const double tol, const size_t max_iter,
			       const vector<double> &vals_hist){
  // EMhist will be the working histogram to do EM
  vector<double> EMhist(vals_hist);
  const double vals_size = accumulate(vals_hist.begin(), vals_hist.end(), 0.0);
  // ensure that alpha is in proper range
  assert((alpha > 0));
  assert(finite(alpha));
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();
  double score = 0.0;
  for(size_t i = 0; i < max_iter; i++){
    EMhist[0] = expected_zeros(vals_size);

    mu = compute_mean(EMhist);

    EMhist[0] = 0.0;
    score = trunc_log_L(EMhist);
    error = fabs((score - prev_score)/score);
    if(error < tol)
      break;
    prev_score = score;
  }

				      return score;
}

double 
ZTNBD::expected_inverse_sum(const double initial_distinct,
                            const double t){
  const double alph = get_alpha();
  const double m = get_mu();
  const double prob_zero = exp(-log(1 + alph*m)/alph);
  return exp(log(initial_distinct) - log(1 - prob_zero) + 
	     log(1 - exp(-log(1 + alpha*m*t)/alph)));
}

double
ZTNBD::expected_distinct(const double initial_distinct,
			 const double t){
  const double prob_zero = exp(-log(1.0 + alpha*mu)/alpha);
  const double prob_zero_t = exp(-log(1.0 + t*alpha*mu)/alpha);
  return exp(log(initial_distinct) + log(1 - prob_zero_t)
	     - log(1 - prob_zero));
}

double
ZTNBD::expected_mincount(const size_t mincount,
			 const double initial_distinct,
			 const double t){
  const double mean_t = get_mu()*t;
  const double alph = get_alpha();
  const double a_inv = 1.0/alph;

  vector<double> probs_mincount;
  for(size_t i = 0; i < mincount; i++)
    probs_mincount.push_back(gsl_sf_lngamma(i + a_inv)
			     - gsl_sf_lngamma(i + 1)
			     - gsl_sf_lngamma(a_inv)
			     + i*(log(alph) + log(mean_t))
			     - (i - a_inv)*log(1.0 + alph*mean_t));

const double lib_size = 
  initial_distinct/(1.0 - exp(-a_inv*log(1.0 + alpha*get_mu())));
 return lib_size*(1 - log_sum_log_vec(probs_mincount, probs_mincount.size()));
}
