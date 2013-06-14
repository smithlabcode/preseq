/*    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith and Timothy Daley
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

#include <numeric>
#include <vector>
#include <string>
#include <iostream>
#include <ostream>
#include <cassert>
#include <limits>
#include <iomanip>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

#include "smithlab_utils.hpp"
#include "newtons_method.hpp"
#include "moment_sequence.hpp"
#include "library_size_estimates.hpp"
#include "ZTNB.hpp"


using std::string;
using std::vector;
using std::max;
using std::cerr;
using std::endl;
using std::setprecision;

using smithlab::log_sum_log_vec;

// Chao (Biometrics 1987) lower bound
double
chao87_lowerbound_unobserved(const vector<double> &counts_hist) {
  assert(counts_hist.size() >= 3);
  return counts_hist[1]*counts_hist[1]/(2.0*counts_hist[2]);
}


// Chao & Lee (JASA 1992) lower bound
double 
cl92_estimate_unobserved(const vector<double> &counts_hist) {
  
  double sample_size = 0.0;
  for(size_t i = 0; i <  counts_hist.size(); i++)
    sample_size += i*counts_hist[i];
  
  const double distinct_vals = 
    accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

  const double estim_coverage = 1.0 - counts_hist[1]/sample_size;
  const double naive_lowerbound = distinct_vals/estim_coverage;
  
  vector<double> log_coeff_var;
  for (size_t i = 2; i < counts_hist.size(); ++i)
    if (counts_hist[i] > 0)
      log_coeff_var.push_back(log(i) + log(i - 1) + log(counts_hist[i]));
  
  const double coeff_var = 
    naive_lowerbound*exp(log_sum_log_vec(log_coeff_var, log_coeff_var.size())
			 -log(sample_size) - log(sample_size - 1)) - 1;
  
  return naive_lowerbound + sample_size*(1.0 - estim_coverage)*coeff_var/estim_coverage;
}

// Chao & Lee (JASA 1992) lower bound
double 
cl92_truncated_estimate_unobserved(const vector<double> &counts_hist,
				   const size_t truncation_count) {
  
  double sample_size = 0.0;
  for(size_t i = 0; i <  counts_hist.size(); i++)
    sample_size += i*counts_hist[i];
  
  const double distinct_vals = 
    accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

  const double estim_coverage = 1.0 - counts_hist[1]/sample_size;
  const double naive_lowerbound = distinct_vals/estim_coverage;
  
  vector<double> log_coeff_var;
  for (size_t i = 2; i < std::min(counts_hist.size(), 
				  truncation_count + 1); ++i)
    if (counts_hist[i] > 0)
      log_coeff_var.push_back(log(i) + log(i - 1) + log(counts_hist[i]));
  
  const double coeff_var = 
    naive_lowerbound*exp(log_sum_log_vec(log_coeff_var, log_coeff_var.size())
			 -log(sample_size) - log(sample_size - 1)) - 1;
  
  return naive_lowerbound + sample_size*(1.0 - estim_coverage)*coeff_var/estim_coverage;
}

///////////////////////////////////////////////
// Methods for Harris lower bounds

//Harris(1959) formula 48
double
harris_3moments_unobserved(const bool VERBOSE,
			   const vector<double> &counts_hist){
  const double m1 = exp(gsl_sf_lnfact(2) + log(counts_hist[2])
			- log(counts_hist[1]));
  const double m2 = exp(gsl_sf_lnfact(3) + log(counts_hist[3])
			- log(counts_hist[1]));
  const double m3 = exp(gsl_sf_lnfact(4) + log(counts_hist[4])
			- log(counts_hist[1]));

  const double sigma2 = m2 - m1*m1;
  const double m3prime = m3 - 3*m1*sigma2 - pow(m1, 3);
  const double y0 = 
    (-m3prime - sqrt(m3prime*m3prime + 4*pow(sigma2, 3)))/(2*sigma2);

  // x's are points, lambda's are weights
  const double x1 = (m1*y0 + sigma2)/y0;
  const double x2 = m1 - y0;
  const double lambda1 = y0*y0/(y0*y0 + sigma2);
  const double lambda2 = 1.0 - lambda1;
  if(VERBOSE){
    cerr << "POINTS = " << setprecision(16) << x2 << ", " << x1 << endl;
    cerr << "WEIGHTS = " << setprecision(16) << lambda2 << ", " << lambda1 << endl;
  }

  return counts_hist[1]*(lambda1/x1 + lambda2/x2);
}




// generate random initial starting position
static void
generate_rand_initial(const size_t dim, const double values_sum,
		      vector<double> &initial_lambdas,
		      vector<double> &initial_xs){
  initial_lambdas.clear();
  initial_xs.clear();

  // random starting x's w/ random weights to start
  for(size_t i = 0; i < dim; i++){
    initial_xs.push_back(static_cast<double>(rand() % static_cast<int>(values_sum)));
    initial_lambdas.push_back(static_cast<double>(rand()));
  }

  // normalize lambdas
  const double lambdas_sum = accumulate(initial_lambdas.begin(), 
					initial_lambdas.end(), 0.0);
  for(size_t i = 0; i < initial_lambdas.size(); i++)
    initial_lambdas[i] = initial_lambdas[i]/lambdas_sum;

  sort(initial_xs.begin(), initial_xs.end());
}

double
harris_newton_unobserved(const bool VERBOSE,
			 const vector<double> &counts_hist,
			 const double tolerance,
			 const size_t max_iter,
			 const size_t n_points){
  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];
  const size_t depth = 2*n_points;
  // depth corresponds to k in Harris
  vector<double> measure_moments;
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < depth; i++)
    measure_moments.push_back(exp(gsl_sf_lnfact(i + 1)
				  + log(counts_hist[i + 1])
				  - log(counts_hist[1])));

  const size_t dim = measure_moments.size()/2;
  vector<double> initial_xs, initial_lambdas;


  vector<double> root_lambdas, root_xs;
  bool CONVERGENCE = false;
  size_t iterations = 0;
  bool CONTINUE = true;
  do{
    if(VERBOSE)
      cerr << "iteration # " << iterations << "\t" << CONVERGENCE << endl;

    generate_rand_initial(dim, values_sum, initial_lambdas, initial_xs);

    if(VERBOSE){
      cerr << "initial xs : ";
      for(size_t i = 0; i < initial_xs.size(); i++)
	cerr << initial_xs[i] << ", ";
      cerr << endl;
      cerr << "initial lambdas : ";
      for(size_t i = 0; i < initial_lambdas.size(); i++)
	cerr << initial_lambdas[i] << ", ";
      cerr << endl;
    }

    // one full iteration of my modified newtons method
    CONVERGENCE = 
      modified_newtons_method(VERBOSE, initial_lambdas, initial_xs, 
			      measure_moments, values_sum, tolerance, 
			      max_iter, root_lambdas, root_xs);
    iterations++;
    // if newtons method did not converge, 
    // start over from random starting position
    CONTINUE = (iterations < max_iter) && !CONVERGENCE;
  }while(CONTINUE);

  double lower_bound = 0.0;
  if(CONVERGENCE){
    if(VERBOSE){
      cerr << "x = ";
      for(size_t i = 0; i < root_xs.size(); i++)
	cerr << root_xs[i] << ", ";
      cerr << endl;
      cerr << "lambdas = ";
      for(size_t i = 0; i < root_lambdas.size(); i++)
	cerr << root_lambdas[i] << ", ";
      cerr << endl;
    }
    for(size_t i = 0; i < root_xs.size(); i++)
      lower_bound += exp(log(counts_hist[1])
			 + log(root_lambdas[i])
			 - log(root_xs[i]));

  }
  else
    cerr << "did not converge" << endl;
  

  return lower_bound;
}


///////////////////////////////////////////////////////
// Quadrature Methods to Estimate Library Size

double
quadrature_unobserved(const bool VERBOSE,
		      const std::vector<double> &counts_hist,
		      const double tol,
		      const size_t max_iter,
		      size_t &n_points){
  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];

  size_t counts_before_first_zero = 1;
  while ((counts_before_first_zero < counts_hist.size()) &&
	  (counts_hist[counts_before_first_zero] > 0))
    ++counts_before_first_zero;
  if(2*n_points > counts_before_first_zero - 2)
    n_points = 
      static_cast<size_t>(floor((counts_before_first_zero - 2)/2));
  
   
  // initialize moments, 0th moment is 1
  vector<double> measure_moments(1, 1.0);
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < 2*n_points; i++)
    measure_moments.push_back(exp(gsl_sf_lnfact(i + 2)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));

  MomentSequence mom_seq(measure_moments);

  vector<double> points, weights;
  mom_seq.QR_quadrature_rules(VERBOSE, n_points,
			      tol,  max_iter, points, weights);

  n_points = points.size();

  const double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
  if(weights_sum != 1){
    cerr << "weights sum = " << weights_sum << endl;
    for(size_t i = 0; i < weights.size(); i++)
      weights[i] = weights[i]/weights_sum;
  }


  // En_1 * int_0^\infty 1/x d \nu (x)
  double estimated_integral = 0.0;
  for(size_t i = 0; i < points.size(); i++)
    estimated_integral += counts_hist[1]*weights[i]/points[i];

  if(VERBOSE){
    cerr << "points = ";
    for(size_t i = 0; i < points.size(); i++)
      cerr << setprecision(16) << points[i] << ", ";
    cerr << endl;

    cerr << "weights = ";
    for(size_t i = 0; i < weights.size(); i++)
      cerr << setprecision(16) << weights[i] << ", ";
    cerr << endl;

    //    cerr << "estimated lib size = " << estimated_integral << endl;
  }


  return estimated_integral;
}


double
quadrature_mean3term_unobserved(const bool VERBOSE,
				const vector<double> &counts_hist,
				const double tol,
				const size_t max_iter,
				const size_t bootstraps,
				size_t &n_points){
  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];

  size_t counts_before_first_zero = 1;
  while ((counts_before_first_zero < counts_hist.size()) &&
	  (counts_hist[counts_before_first_zero] > 0))
    ++counts_before_first_zero;
  if(2*n_points > counts_before_first_zero - 2)
    n_points = 
      static_cast<size_t>(floor((counts_before_first_zero - 2)/2));
  
   
  // initialize moments, 0th moment is 1
  vector<double> measure_moments(1, 1.0);
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < 2*n_points; i++)
    measure_moments.push_back(exp(gsl_sf_lnfact(i + 2)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));

  MomentSequence mom_seq(measure_moments);

  mom_seq.set_mean_3term_recurrence(VERBOSE, counts_hist, 
				    bootstraps, n_points);

  vector<double> points, weights;
  mom_seq.QR_quadrature_rules(VERBOSE, n_points,
			      tol,  max_iter, points, weights);

  n_points = points.size();

  const double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
  if(weights_sum != 1){
    cerr << "weights sum = " << weights_sum << endl;
    for(size_t i = 0; i < weights.size(); i++)
      weights[i] = weights[i]/weights_sum;
  }


  // En_1 * int_0^\infty 1/x d \nu (x)
  double estimated_integral = 0.0;
  for(size_t i = 0; i < points.size(); i++)
    estimated_integral += counts_hist[1]*weights[i]/points[i];

  if(VERBOSE){
    cerr << "points = ";
    for(size_t i = 0; i < points.size(); i++)
      cerr << setprecision(16) << points[i] << ", ";
    cerr << endl;

    cerr << "weights = ";
    for(size_t i = 0; i < weights.size(); i++)
      cerr << setprecision(16) << weights[i] << ", ";
    cerr << endl;

    //    cerr << "estimated lib size = " << estimated_integral << endl;
  }

  return estimated_integral;
}


// need a better name
// computes quadrature rules using points determined 
// assuming NegBin count distribution 
// then computes weights to match moment conditions
double
NegBin_quadrature_unobserved(const bool VERBOSE,
			     const std::vector<double> &counts_hist,
			     const double tol,
			     const size_t max_iter,
			     size_t &n_points){
  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];

  size_t counts_before_first_zero = 1;
  while ((counts_before_first_zero < counts_hist.size()) &&
	  (counts_hist[counts_before_first_zero] > 0))
    ++counts_before_first_zero;
  if(2*n_points > counts_before_first_zero - 2)
    n_points = 
      static_cast<size_t>(floor((counts_before_first_zero - 2)/2));
  
   
  // initialize moments, 0th moment is 1
  vector<double> measure_moments(1, 1.0);
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < 2*n_points; i++)
    measure_moments.push_back(exp(gsl_sf_lnfact(i + 2)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));

  MomentSequence mom_seq(measure_moments);
  vector<double> NBhist(counts_hist);
  ZTNBD distro(1.0, 1.0);
  distro.EM_estim_params(tol, max_iter, NBhist);
  if(VERBOSE){
    cerr << "estimated mu    = " << distro.get_mu() << endl;
    cerr << "estimated alpha = " << distro.get_alpha() << endl;
  }

  vector<double> points, weights;
  mom_seq.NegBin_quadrature_rules(VERBOSE, n_points,
				  tol,  max_iter, 
				  distro.get_mu(),
				  distro.get_alpha(),
				  points, weights);

  n_points = points.size();

  const double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
  if(weights_sum != 1){
    cerr << "weights sum = " << weights_sum << endl;
    for(size_t i = 0; i < weights.size(); i++)
      weights[i] = weights[i]/weights_sum;
  }


  // En_1 * int_0^\infty 1/x d \nu (x)
  double estimated_integral = 0.0;
  for(size_t i = 0; i < points.size(); i++)
    estimated_integral += counts_hist[1]*weights[i]/points[i];

  if(VERBOSE){
    cerr << "points = ";
    for(size_t i = 0; i < points.size(); i++)
      cerr << setprecision(16) << points[i] << ", ";
    cerr << endl;

    cerr << "weights = ";
    for(size_t i = 0; i < weights.size(); i++)
      cerr << setprecision(16) << weights[i] << ", ";
    cerr << endl;

    //    cerr << "estimated lib size = " << estimated_integral << endl;
  }


  return estimated_integral;
}
