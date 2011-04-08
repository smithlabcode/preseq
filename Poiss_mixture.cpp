/*    Poiss_mixture:
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


static inline double 
compute_mean(const vector<size_t> &vals_hist){
  double mean = 0.0;
  const double vals_size = static_cast<double>(
                       accumulate(vals_hist.begin(), 
                                  vals_hist.end(), 0));
  for(size_t i = 0; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}


static inline double
compute_weighted_mean_hist(const vector<size_t> &vals_hist,
			   const vector<double> &probs){
  double mean = 0.0;
  double probs_weight = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++){
    mean += i*probs[i]*vals_hist[i];
    probs_weight += vals_hist[i]*probs[i];
  }
  return(mean/probs_weight);
}

void 
Poiss::estim_param(const vector<size_t> &vals_hist){
  lambda = compute_mean(vals_hist);
  assert(finite(lambda));
}

void 
Poiss::estim_param(const vector<size_t> &vals_hist, 
                   const vector<double> &probs){
  lambda = compute_weighted_mean_hist(vals_hist,probs);

}


static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

static inline double
lambda_score_funct(const double mean, const double lambda){
  return(lambda + mean*exp(-lambda) - mean);
}

void
ZTP::trunc_estim_param_bisec(const double mean, const double tol){

  double lambda_low = mean-1;
  double lambda_high = mean;
  double lambda_mid = mean - 0.5;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while(movement(lambda_high, lambda_low) > tol &&
         diff > tol){
    
    lambda_mid = (lambda_low + lambda_high)/2;
    
    const double mid_val = lambda_score_funct(mean, lambda_mid);
    
    if (mid_val < 0) lambda_low = lambda_mid;
    else lambda_high = lambda_mid;
    
    diff = fabs((prev_val - mid_val)/max(mid_val, prev_val));
    
    prev_val = mid_val;
  }
  set_lambda(lambda_mid);
}

double 
ZTP::expected_inverse_sum(const size_t sample_size,
			  const size_t sum){
    const double lamb = get_lambda();
    const double num = exp(log(sum)-log(sample_size)+log(1-exp(-lamb)));
    return(exp(log(sample_size)-log(1-exp(-lamb))
	       +log(1-exp(-num))));
}



double
ZTP_mixture::trunc_expectation_step(const vector<size_t> &vals_hist,
                         vector< vector<double> > &probs){
  double score = 0.0;

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

    score += vals_hist[i]*log_denom;
  }
  return(score);
}

void
ZTP_mixture::trunc_calculate_mixing(const vector<size_t> &vals_hist,
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
ZTP_mixture::trunc_max_step(const vector<size_t> &vals_hist,
			    const vector< vector<double> > &probs, 
			    const double tol){

  for(size_t i = 0; i < distros.size(); i++){
    const double mean = compute_weighted_mean_hist(vals_hist, probs[i]);
    distros[i].trunc_estim_param_bisec(mean, tol);
  }

  trunc_calculate_mixing(vals_hist, probs);

}
// recalculate parameters, lambda and mixing_j = average(probs[j]) //

double 
ZTP_mixture::trunc_log_L(const vector<size_t> &vals_hist){
  double logL = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++){
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

double
ZTP_mixture::Fisher_obs_info_2nd_theta(const vector< vector<double> > &probs,
				       const vector<size_t> &vals_hist,
				       const size_t indx){
  /* I(complete) = sum_j 1(x_j \in indx)/mixing[i]^2
                         + 1(x_j \in last)/mixing.back()^2
  I(X| complete) = 
    sum_j (1(x_j \in i)/mixing[i] - 1(x_j \in last)/mixing.back)^2
    - (probs[indx][j]/mixing[i] - probs[last][j]/mixing.back())^2
  I(X) = I(complete) - I(X | complete)
  */
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_prob = log(probs[indx][i]);
      const double log_last_prob = log(probs[probs.size()-1][i]);
      if(finite(log_prob)){
	pos_log_vec.push_back(log_vals_hist + 2*log_prob
			      -2*log(mixing[indx]));
      }
      if(finite(log_last_prob)){
	pos_log_vec.push_back(log_vals_hist
			      +2*log_last_prob
			      -2*log(mixing.back()));
      }
      if(finite(log_prob) &&
	 finite(log_last_prob)){
	neg_log_vec.push_back(log_vals_hist+log(2)+log_prob
			      +log_last_prob
			      -log(mixing[indx])
			      -log(mixing.back()));
      }
    }
  }
   
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double
ZTP_mixture::Fisher_obs_info_mixed_theta(
                                 const vector< vector<double> > &probs,
				 const vector<size_t> &vals_hist,
				 const size_t indx1,
				 const size_t indx2){
  /* I(complete) = \sum_j probs[last]/mixing.back()^2
     I(X | complete) = \sum_j probs[last]/mixing.back()^2
        -probs[indx1][j]*probs[indx2][j]/(mixing[indx1]*mixing[indx2])
        +(probs[last][j]/mixing.back())*(probs[indx1][j]/mixing[indx1]
                                         +probs[indx2][j]/mixing[indx2])
        -(probs[last][j]/mixing.back())^2
  */
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
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
ZTP_mixture::Fisher_obs_info_2nd_lambda(
                                const vector< vector<double> > &probs,
				const vector<size_t> &vals_hist,
				const size_t indx){
  /* I(complete) = \sum_j 1(x_j \in i)*x_j/lambda_i^2
                          +1(x_j \in i)*exp(lambda_i)/(exp(lambda_i)-1)
                          -1(x_j \in i)*exp(2*lambda_i)/(exp(lambda_i)-1)^2
     I(X|complete) = 
      \sum_j probs[i][j]*(x_j/lambda_i -(e^lambda_i)/(e^lambda_i -1))^2
             -probs[i][j]^2*(x_j/lambda_i -(e^lambda_i)/(e^lambda_i -1))^2
  */
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  const double lambda = distros[indx].get_lambda();
  for(size_t i =1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_prob = log(probs[indx][i]);
      if(finite(log_prob)){
	pos_log_vec.push_back(log_vals_hist 
			      +log_prob + log(i)
			      -2*log(lambda));
	pos_log_vec.push_back(log_vals_hist
			      +log_prob
			      -log(1-exp(-lambda)));
	pos_log_vec.push_back(log_vals_hist
			      +log_prob 
			      +2*log(fabs(i/lambda - 1/(1-exp(-lambda)))));
        neg_log_vec.push_back(log_vals_hist
			      +log_prob
			      -2*log(1-exp(-lambda)));
	neg_log_vec.push_back(log_vals_hist
			      +2*log_prob 
			      +2*log(fabs(i/lambda - 1/(1-exp(-lambda)))));
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double
ZTP_mixture::Fisher_obs_info_mixed_lambda(
                                  const vector< vector<double> > &probs,
				  const vector<size_t> &vals_hist,
				  const size_t indx1,
				  const size_t indx2){
  /* I(complete) = 0 */
  const double lambda1 = distros[indx1].get_lambda();
  const double lambda2 = distros[indx2].get_lambda();
  vector<double> summands;
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      double deriv_1 = exp(log(i)-log(lambda1)) - 1/(1-exp(-lambda1));
      double deriv_2 = exp(log(i)-log(lambda2)) - 1/(1-exp(-lambda2));
      summands.push_back(vals_hist[i]*probs[indx1][i]*probs[indx2][i]*
			 deriv_1*deriv_2);
    }
  }
  return(accumulate(summands.begin(), summands.end(), 0.0));
}

double
ZTP_mixture::Fisher_obs_info_mixed_same_indx( 
                                      const vector< vector<double> > &probs,
				      const vector<size_t> &vals_hist,
				      const size_t indx){
  /*I(complete) = 0 */
  const double lambda = distros[indx].get_lambda();
  const double theta = mixing[indx];
  const double last_mixing = mixing.back();
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  for(size_t i =1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_prob = log(probs[indx][i]);
      const double log_last_prob = log(probs[probs.size()-1][i]);
      if(finite(log_prob)){
	pos_log_vec.push_back(log_vals_hist+log_prob-log(theta)
			      -log(1-exp(-lambda)));
	pos_log_vec.push_back(log_vals_hist+2*log_prob+log(i)
			      -log(lambda)-log(theta));
	neg_log_vec.push_back(log_vals_hist+log_prob+log(i)
			      -log(lambda)-log(theta));
	neg_log_vec.push_back(log_vals_hist+2*log_prob-log(theta)
			      -log(1-exp(-lambda)));
	if(finite(log_last_prob)){
	  pos_log_vec.push_back(log_vals_hist+log_prob+log_last_prob
				-log(last_mixing)-log(1-exp(-lambda)));
	  neg_log_vec.push_back(log_vals_hist+log_prob+log_last_prob
				+log(i)-log(lambda)-log(last_mixing));
	}
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}


double 
ZTP_mixture::Fisher_obs_info_mixed_mixed_indx(
                                      const vector< vector<double> > &probs,
				      const vector<size_t> &vals_hist,
				      const size_t lambda_indx,
				      const size_t theta_indx){
  /*I(complete) = 0 */
  const double theta = mixing[theta_indx];
  const double lambda = distros[lambda_indx].get_lambda();
  const double last_mixing = mixing.back();
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_prob_lambda = log(probs[lambda_indx][i]);
      const double log_prob_theta = log(probs[theta_indx][i]);
      const double log_last_prob = log(probs[probs.size()-1][i]);
      if(finite(log_prob_lambda) && finite(log_prob_theta)){
	pos_log_vec.push_back(log_vals_hist + log_prob_lambda 
			      + log_prob_theta +log(i) - log(lambda)
			      -log(theta));
	neg_log_vec.push_back(log_vals_hist + log_prob_lambda
			      +log_prob_theta - log(theta)
			      -log(1-exp(-lambda)));
      }
      if(finite(log_prob_lambda) && finite(log_last_prob)){
	pos_log_vec.push_back(log_vals_hist + log_prob_lambda
			      +log_last_prob -log(last_mixing)
			      -log(1-exp(-lambda)));
	neg_log_vec.push_back(log_vals_hist+log_prob_lambda
			      -log_last_prob+log(i)-log(lambda)
			      -log(last_mixing));
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}

double 
ZTP_mixture::Fisher_obs_info_mixed_last_indx(
                                     const vector< vector<double> > &probs,
				     const vector<size_t> &vals_hist,
				     const size_t theta_indx){
  /* indx of lambda  = K
     I(complete) = 0 */
  const double theta = mixing[theta_indx];
  const double lambda = distros[distros.size()-1].get_lambda();
  const double last_mixing = mixing.back();
  vector<double> pos_log_vec;
  vector<double> neg_log_vec;
  for(size_t i = 1; i < vals_hist.size(); i++){
    if(vals_hist[i] > 0){
      const double log_vals_hist = log(vals_hist[i]);
      const double log_prob_theta = log(probs[theta_indx][i]);
      const double log_last_prob = log(probs[probs.size()-1][i]);  
      if(finite(log_last_prob)){
	neg_log_vec.push_back(log_vals_hist + log_last_prob
			      +log(i) - log(lambda) - log(last_mixing));
	pos_log_vec.push_back(log_vals_hist + log_last_prob
			      -log(last_mixing) - log(1-exp(-lambda)));
	neg_log_vec.push_back(log_vals_hist + 2*log_last_prob
			      -log(last_mixing) + log(i) -  log(lambda));
	pos_log_vec.push_back(log_vals_hist + 2*log_last_prob
			      - log(last_mixing) - log(1-exp(-lambda)));
	if(finite(log_prob_theta)){
	  pos_log_vec.push_back(log_vals_hist + log_last_prob
				+log_prob_theta + log(i)
				-log(lambda) - log(theta));
	  neg_log_vec.push_back(log_vals_hist + log_last_prob
				+log_prob_theta -log(theta)
				-log(1-exp(-lambda)));
	}
      }
    }
  }
  return(exp(log_sum_log_vec(pos_log_vec, pos_log_vec.size()))
	 -exp(log_sum_log_vec(neg_log_vec, neg_log_vec.size())));
}


void
ZTP_mixture::compute_observed_Fisher_info(const vector<size_t> &vals_hist,
					    const vector< vector<double> >
					    &probs){
  const size_t number_states = distros.size();
  for(size_t i = 0; i < number_states-1; i++){
    for(size_t j = 0; j < number_states-1; j++){
      if(i != j){
	Fisher_info[i][j] = 
	  Fisher_obs_info_mixed_theta(probs, vals_hist, i,j);
      }
      else{
	Fisher_info[i][j] = 
	  Fisher_obs_info_2nd_theta(probs, vals_hist, i);
      }
    }
  }
  for(size_t i = number_states-1; i < 2*number_states-1; i++){
    for(size_t j = number_states-1; j < 2*number_states-1; j++){
      if(i != j){
	Fisher_info[i][j] = 
	  Fisher_obs_info_mixed_lambda(probs,
					    vals_hist, i-number_states+1,
					    j-number_states+1);
      }
      else{
	Fisher_info[i][j] = 
	  Fisher_obs_info_2nd_lambda(probs, vals_hist, i-number_states+1);
      }
    }
  }
  for(size_t i = 0;  i < number_states-1; i++){
    for(size_t j = 0; j < number_states-1; j++){
      if(i != j){
	Fisher_info[i][j+number_states-1] = 
	  Fisher_obs_info_mixed_mixed_indx(probs, vals_hist, j, i);
	Fisher_info[j+number_states-1][i] = 
	  Fisher_info[i][j+number_states-1];
      }
      else{
	Fisher_info[i][j+number_states-1] = 
	  Fisher_obs_info_mixed_same_indx(probs, vals_hist, i);
	Fisher_info[j+number_states-1][i] =
	  Fisher_info[i][j+number_states-1];
      }
    }
  }
  for(size_t i=0; i< number_states-1; i++){
    Fisher_info[i][2*number_states-2] = 
      Fisher_obs_info_mixed_last_indx(probs, vals_hist, i);
    Fisher_info[2*number_states-2][i] = 
      Fisher_info[i][2*number_states-2];
  }

}

double 
ZTP_mixture::EM_mix_resolve(const vector<size_t> &vals_hist,
			    const double tol,const size_t max_iter){
  
  const size_t number_states = distros.size();
  double probs_starting_val = 1/static_cast<double>(number_states);
  vector< vector<double> > probs(number_states, 
                                 vector<double>(vals_hist.size(),
                                                probs_starting_val));

  double score = 0.0;
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_iter; ++i){

    score = trunc_expectation_step(vals_hist, probs);
    trunc_max_step(vals_hist, probs, tol);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;


  }
  /* compute Fisher info */
  compute_observed_Fisher_info(vals_hist, probs);


  return(trunc_log_L(vals_hist));
}

double
Poiss_mixture::expectation_step(const vector<size_t> &vals_hist,
				vector< vector<double> > &probs){
  double score = 0.0;

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

    score += vals_hist[i]*log_denom;
  }
  return(score);
}

void
Poiss_mixture::calculate_mixing(const vector<size_t> &vals_hist,
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
Poiss_mixture::maximization_step(const vector<size_t> &vals_hist,
			    const vector< vector<double> > &probs){

  for(size_t i = 0; i < distros.size(); i++){
    distros[i].estim_param(vals_hist);
  }

  calculate_mixing(vals_hist, probs);

}
// recalculate parameters, lambda and mixing_j = average(probs[j]) //

double 
Poiss_mixture::log_L(const vector<size_t> &vals_hist){
  double logL = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++){
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

double 
Poiss_mixture::EM_mix_resolve(const vector<size_t> &vals_hist,
			    const double tol,const size_t max_iter){
  
  const size_t number_states = distros.size();
  double probs_starting_val = 1/static_cast<double>(number_states);
  vector< vector<double> > probs(number_states, 
                                 vector<double>(vals_hist.size(),
                                                probs_starting_val));

  double score = 0.0;
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_iter; ++i){

    score = expectation_step(vals_hist, probs);
    maximization_step(vals_hist, probs);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;


  }

  return(log_L(vals_hist));
}

double
ZTP_mixture::expected_inverse_sum(const size_t sample_size,
				  const size_t sum){
  const vector<ZTP> dist = get_distros();
  const vector<double> mix = get_mixing();
  vector<double> lambdas;
  for(size_t i = 0; i < dist.size(); i++)
    lambdas.push_back(dist[i].get_lambda());

  vector<double> log_mean_vec;
  for(size_t i = 0; i < lambdas.size(); i++){
    log_mean_vec.push_back(log(mix[i]) + log(lambdas[i])
			   -log(1-exp(-lambdas[i])));
  }
  const double mean = exp(log_sum_log_vec(log_mean_vec,
					  log_mean_vec.size()));
  vector<double> log_terms_vec;
  for(size_t i = 0; i < lambdas.size(); i++){
    const double num = 
      exp(log(sum) - log(sample_size) +log(lambdas[i]) - log(mean));
    log_terms_vec.push_back(log(sample_size) + log(mix[i]) 
			    -log(1-exp(-lambdas[i]))+log(1-exp(-num)));
  }
  return(exp(log_sum_log_vec(log_terms_vec, log_terms_vec.size())));
}
