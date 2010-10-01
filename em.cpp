#include "em.hpp"

#include <memory>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <limits>

#include <gsl/gsl_randist.h>

using std::vector;
using std::cerr;
using std::endl;

inline double
log_sum_log_vec(const vector<double> &vals, size_t limit) {
  const vector<double>::const_iterator x = 
    max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
      assert(finite(sum));
    }
  }
  return max_val + log(sum);
}

inline static double
log_sum_log(const double p, const double q) {
  if (p == 0)
    return q;
  else if (q == 0)
    return p;
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


inline static double
log_sum_log(const double p, const double q, const double r) {
  if (p == 0) return log_sum_log(q, r);
  else if (q == 0) return log_sum_log(p, r);
  else if (r == 0) return log_sum_log(p, r);
  else return log_sum_log(p, log_sum_log(q, r));
}

static double
expectation_step(const vector<double> &values, 
		 const vector<double> &mixing,
		 const double &first_lambda, 
		 const double &second_lambda, 
		 const double &third_lambda,
		 vector<double> &first_probs, 
		 vector<double> &second_probs, 
		 vector<double> &third_probs) {
  
/*  vector<double>::const_iterator x_idx = values.begin();
  const vector<double>::const_iterator x_lim = values.end();*/
  
  double score = 0;
  
  const double first_log_mixing = log(mixing[0]);
  assert(finite(first_log_mixing));
  const double second_log_mixing = log(mixing[1]);
  assert(finite(second_log_mixing));
  const double third_log_mixing = log(mixing[2]);
  assert(finite(third_log_mixing));
  
  for (size_t i = 0; i < values.size(); ++i) {
    
    const double first_part = first_log_mixing + log(gsl_ran_poisson_pdf(values[i], first_lambda));
    assert(finite(first_part));
    
    const double second_part = second_log_mixing + log(gsl_ran_poisson_pdf(values[i], second_lambda));
    assert(finite(second_part));

    const double third_part = third_log_mixing + log(gsl_ran_poisson_pdf(values[i], third_lambda));
    assert(finite(first_part));
    
    const double denom = log_sum_log(first_part, second_part, third_part);
    assert(finite(denom));
    
    first_probs[i] = exp(first_part - denom);
    second_probs[i] = exp(second_part - denom);
    third_probs[i] = exp(third_part - denom);
    
    score += denom;
  }
  
  return score;
}


static inline double
weighted_mean(const vector<double> &weights, 
	      const vector<double> &values) {
  return std::inner_product(weights.begin(), weights.end(), 
			    values.begin(), 0.0)/
    std::accumulate(weights.begin(), weights.end(), 0.0);
}

static void
maximization_step(const vector<double> &values, 
		  const vector<double> &first_probs,
		  const vector<double> &second_probs,
		  const vector<double> &third_probs,
		  vector<double> &mixing,
		  double &first_lambda, 
		  double &second_lambda, 
		  double &third_lambda) {

  first_lambda = weighted_mean(first_probs, values);
  second_lambda = weighted_mean(second_probs, values);
  third_lambda = weighted_mean(third_probs, values);
  
  vector<double> log_first_probs(first_probs);
  vector<double> log_second_probs(second_probs);
  vector<double> log_third_probs(third_probs);
  for (size_t i = 0; i < log_first_probs.size(); ++i) {
    log_first_probs[i] = log(log_first_probs[i]);
    log_second_probs[i] = log(log_second_probs[i]);
    log_third_probs[i] = log(log_third_probs[i]);
  }    
  
  mixing[0] = log_sum_log_vec(log_first_probs, log_first_probs.size());
  mixing[1] = log_sum_log_vec(log_second_probs, log_second_probs.size());
  mixing[2] = log_sum_log_vec(log_third_probs, log_third_probs.size());
  const double mix_sum = log_sum_log(mixing[0], mixing[1], mixing[2]);
  
  mixing[0] = exp(mixing[0] - mix_sum);
  mixing[1] = exp(mixing[1] - mix_sum);
  mixing[2] = exp(mixing[2] - mix_sum);
}


void
EM(const vector<double> &values,
   const size_t max_iterations, const double tolerance, 
   const bool VERBOSE,
   double &first_lambda, double &second_lambda, double &third_lambda,
   vector<double> &mixing) {
  
  vector<double> first_probs(values.size(), 0);
  vector<double> second_probs(values.size(), 0);
  vector<double> third_probs(values.size(), 0);
  
  mixing = vector<double>(3, 1.0/3.0);
  
  if (VERBOSE)
    cerr << endl << std::setw(10) << "DELTA"
	 << std::setw(14) << "(PARAMS,MIX)" << endl;
  
  // Do the expectation maximization
  double prev_score = std::numeric_limits<double>::max();
  for (size_t itr = 0; itr < max_iterations; ++itr) {
    const double score = 
      expectation_step(values, mixing, first_lambda, second_lambda, third_lambda, 
		       first_probs, second_probs, third_probs);
    maximization_step(values, first_probs, second_probs, third_probs, 
		      mixing, first_lambda, second_lambda, third_lambda);
    if (VERBOSE) {
      cerr << std::setw(10) << std::setprecision(4) 
	   << (prev_score - score)/prev_score << endl
	   << std::setw(14) << first_lambda << " " 
	   << std::setw(10) << mixing[0] << endl
	   << std::setw(14) << second_lambda << " " 
	   << std::setw(10) << mixing[1] << endl
	   << std::setw(14) << third_lambda << " " 
	   << std::setw(10) << mixing[2] << endl
	   << endl;
    }
    if (std::fabs((prev_score - score)/prev_score) < tolerance)
      break;
    prev_score = score;
  }
}
