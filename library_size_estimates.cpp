/*    Copyright (C) 2012 University of Southern California and
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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_gamma.h>

#include "smithlab_utils.hpp"

#include "library_size_estimates.hpp"
#include "newtons_method.hpp"
#include "gaussian_quadrature.hpp"

using std::string;
using std::vector;
using std::max;
using std::cerr;
using std::endl;

using smithlab::log_sum_log_vec;

static const size_t MIN_ALLOWED_DEGREE = 6;

/* compute the upperbound on library size by noting that p(x)/q(x) is
 * a constant in the limit if p & q are of equal degree. furthermore
 * if p & q are equal degree then the library_yield estimates are
 * liberal, i.e. they overestimate the library_yield on average
 */
/*
double
upperbound_librarysize(const bool VERBOSE, const vector<double> &counts_hist, 
		       size_t max_terms) {
  // need max_terms = L + M + 1 to be even so that L + M is odd so
  // that we can take lim_{t \to \infty} [L+1, M]
  
  if (max_terms % 2 == 1)
    --max_terms;
  
  vector<double> ps_coeffs;
  for (size_t i = 0; i < max_terms; ++i)
    ps_coeffs.push_back(pow(-1, i + 2)*counts_hist[i + 1]);
  
  for (; max_terms >= MIN_ALLOWED_DEGREE; max_terms -= 2) {
    const size_t numer_degree = max_terms/2;  
    const size_t denom_degree = max_terms - numer_degree;
    assert(numer_degree == denom_degree);
    // numer_degree = degree(p)+1, denom_degree = degree(q); and
    // library_yield = x*p(x)/q(x), ensure that degree(x*p)=degree(q)
    
    // consider upper bound if pade approx is acceptable
    vector<double> numers, denoms;
    const bool accept_approx = 
      compute_pade_coeffs(ps_coeffs, numer_degree, denom_degree, numers, denoms); 
    
    // lim(xp(x)/q(x)) = p_{numer_size-1}/q_{denom_size} coefficients
    // are in order of degree
    const double upper_bound = numers.back()/denoms.back();
    if (accept_approx && upper_bound > 0.0 && std::isfinite(upper_bound)){
      if(VERBOSE){
	cerr << "UPPER_BOUND_NUMERATOR_COEFFICIENTS" <<  endl;
	for(size_t j = 0; j < numers.size(); j++)
	  cerr << numers[j] << endl;
	cerr << "UPPER_BOUND_DENOMINATOR_COEFFICIENTS" << endl;
	for(size_t j = 0; j < denoms.size(); j++)
	  cerr << denoms[j] << endl;
      }

      // use highest order acceptable approximation
      return upper_bound;

    }
  }
  return -std::numeric_limits<double>::max();
}
*/

// Chao (Biometrics 1987) lower bound
double
chao87_lowerbound_librarysize(const vector<double> &counts_hist) {
  assert(counts_hist.size() >= 2);
  return accumulate(counts_hist.begin(), counts_hist.end(), 0.0) +
    counts_hist[1]*counts_hist[1]/(2.0*counts_hist[2]);
}


// Chao & Lee (JASA 1992) lower bound
double 
cl92_estimate_librarysize(const vector<double> &counts_hist) {
  
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
  
  return naive_lowerbound + counts_hist[1]*coeff_var/estim_coverage;
}

//////////////////////////////////////////////////
// Harris (1959) lower bound
/*
struct pars{
  // observed moments
  vector<double> moments;
  // number of captures, i.e. mapped reads
  double N;
};

static void
set_lambdas(const gsl_vector *lambdas_xs,
	    const size_t dim,
	    vector<double> &lambdas){
  lambdas.clear();
  // lambdas are first dim + 1 terms
  for(size_t i = 0; i < dim; i++)
    lambdas.push_back(gsl_vector_get(lambdas_xs, i));
  lambdas.push_back(gsl_vector_get(lambdas_xs, dim));
}

static void
set_xs(const gsl_vector *lambdas_xs,
       const size_t dim,
       vector<double> &xs){
  xs.clear();
  // x's are last dim terms
  for(size_t  i = dim + 1; i < lambdas_xs->size; i++)
    xs.push_back(gsl_vector_get(lambdas_xs, i));
}

// lambdas_xs is lambdas followed by x's
// last x = N
// eqns: sum_i lambda_i = 1
// sum_i lambda_i x_i = mu_1
// sum_i lambda_i x_i^2 = mu_2
// ....
// sum_i lambda_i x_i^k = mu_k
static int
set_system_eqns(const gsl_vector *lambdas_xs, 
		void *param,
		gsl_vector *f){
  const vector<double> mmnts = ((struct pars*) param)->moments;
  const double vals_sum = ((struct pars*) param)->N;

  const size_t dim = mmnts.size()/2;
  assert(lambdas_xs->size == mmnts.size() + 1);
  vector<double> lambdas, xs;

  set_lambdas(lambdas_xs, dim, lambdas);
  set_xs(lambdas_xs, dim, xs);

  // sum lambda_{i} = 1
  gsl_vector_set(f, 0, accumulate(lambdas.begin(), 
				  lambdas.end(), 0.0) - 1.0);

  for(size_t i = 0; i < mmnts.size(); i++){
    // sum_j lambda_{j}*x_{j}^k = \mu_{k} (last x = N)
    double function_value = - mmnts[i];
    for(size_t j = 0; j < xs.size(); j++)
      function_value += lambdas[j]*pow(xs[j], i + 1);
    function_value += lambdas.back()*pow(vals_sum, i + 1);
    gsl_vector_set(f, i + 1, function_value);
  }

  return GSL_SUCCESS;
}


// compute and set the jacobian
static int
set_system_jacobian(const gsl_vector *lambdas_xs,
		    void *param,
		    gsl_matrix *Jacob){
  const vector<double> mmnts = ((struct pars*) param)->moments;
  const double vals_sum = ((struct pars*) param)->N;
  const size_t dim = mmnts.size()/2;

  vector<double> lambdas, xs;

  set_lambdas(lambdas_xs, dim, lambdas);
  set_xs(lambdas_xs, dim, xs);

  for(size_t j = 0; j < dim + 1; j++)
    gsl_matrix_set(Jacob, 0, j, 1.0);

  // sum lambda_j*x_j^i - mu_k
  for(size_t i = 0; i < 2*dim; i++)
    for(size_t j = 0; j < dim; j++)
      gsl_matrix_set(Jacob, i + 1, j, pow(xs[j], i + 1));

  // last lambda
  for(size_t i = 0; i < 2*dim; i++)
    gsl_matrix_set(Jacob, i + 1, dim, pow(vals_sum, i + 1));

  // now for the x's
  for(size_t i = 0; i < 2*dim; i++)
    for(size_t j = 0; j < dim; j++)
      gsl_matrix_set(Jacob, i + 1, j + dim + 1, 
		     lambdas[j]*(i + 1)*pow(xs[j], i));

  return GSL_SUCCESS;
}


// set function value and jacobian
static int
set_fdf(const gsl_vector *lambdas_xs,
	void *param,
	gsl_vector *f,
	gsl_matrix *Jacob){

  set_system_eqns(lambdas_xs, param, f);
  set_system_jacobian(lambdas_xs, param, Jacob);

  return GSL_SUCCESS;
}

static inline void
print_state(const size_t iter, 
	    const size_t total_dim,
	    const gsl_multiroot_fdfsolver *s){
  cerr << "iter = " << iter << endl;
  double residual = 0.0;
  for(size_t i = 0; i < total_dim; i++)
    residual += pow(gsl_vector_get(s->f, i), 2);
  cerr << "residual = " << residual << endl;
  cerr << "system of equations : ";
  for(size_t i = 0; i < total_dim; i++)
    cerr << gsl_vector_get(s->f, i) << "\t";
  cerr << endl;
}


//Harris (AnnalsMathStat 1959) lower bound, 
// Chao's(1987) bound is a simple approx
// uses gsl's newton's method
double 
harris_lowerbound_librarysize(const vector<double> &counts_hist,
			      const double tolerance,
			      const size_t max_iter,
			      const size_t depth){
  double vals_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    vals_sum += i*counts_hist[i];

  // depth corresponds to k in Harris
  vector<double> measure_moments;
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < depth; i++)
    measure_moments.push_back(exp(gsl_sf_lngamma(i + 3)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));


  // need to use an even number of coefficients
  // odd number gives upper bound ( = infty)
  if(depth % 2 == 1)
    measure_moments.pop_back();

  const gsl_multiroot_fdfsolver_type *solver_type;
  gsl_multiroot_fdfsolver *solver;
  int status = GSL_CONTINUE;
  size_t iter = 0;

  const size_t total_dim = measure_moments.size() + 1;
  const size_t dim = measure_moments.size()/2;
  gsl_vector *lambdas_xs = gsl_vector_alloc(total_dim);
  

  // set initial lambdas to be uniform
  for(size_t i = 0; i < dim + 1; i++)
    gsl_vector_set(lambdas_xs, i, 1.0/(dim + 1.0));

  // set initial x's to be random
  srand(time(NULL));
  vector<double> initial_xs;
  for(size_t i = 0; i < dim; i++)
    initial_xs.push_back(static_cast<double>(rand() % static_cast<int>(vals_sum)));

  sort(initial_xs.begin(), initial_xs.end());
  cerr << "initial x's = ";
  for(size_t i = 0; i < initial_xs.size(); i++)
    cerr << initial_xs[i] << ", ";
  cerr << endl;
  for(size_t i = 0 ; i < dim; i++)
    gsl_vector_set(lambdas_xs, i + dim + 1, initial_xs[i]);


  struct pars p = {measure_moments, vals_sum};

  gsl_multiroot_function_fdf fdf = {&set_system_eqns,
				    &set_system_jacobian,
				    &set_fdf,
				    total_dim, &p};

  solver_type = gsl_multiroot_fdfsolver_hybridsj;
  solver = gsl_multiroot_fdfsolver_alloc(solver_type, total_dim);
  gsl_multiroot_fdfsolver_set(solver, &fdf, lambdas_xs);

  gsl_vector *current_lambdas_xs = gsl_vector_alloc(total_dim);

  // print_state(iter, total_dim, solver);

  while(iter < max_iter & status == GSL_CONTINUE){
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(solver);
    //  print_state(iter, total_dim, solver);
    
    status = gsl_multiroot_test_residual(solver->f, tolerance);
    current_lambdas_xs = gsl_multiroot_fdfsolver_root(solver);
  }

  double lower_bound = 0.0;
  if(status != GSL_CONTINUE){
    vector<double> lambdas;
    vector<double> xs;
    set_lambdas(current_lambdas_xs, measure_moments.size()/2, lambdas);
    set_xs(current_lambdas_xs, measure_moments.size()/2, xs);

    cerr << "lambdas = ";
    for(size_t i = 0; i < lambdas.size(); i++)
      cerr << lambdas[i] << ", ";
    cerr << endl;
    cerr << "x's = ";
    for(size_t i = 0; i < xs.size(); i++)
      cerr << xs[i] << ", ";
    cerr << vals_sum << endl;

    for(size_t i = 0; i < xs.size(); i++)
      lower_bound += counts_hist[1]*lambdas[i]/xs[i];

    lower_bound += counts_hist[1]*lambdas.back()/vals_sum;
  }

  return lower_bound;
}
*/


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
  initial_lambdas.push_back(rand());

  // normalize lambdas
  const double lambdas_sum = accumulate(initial_lambdas.begin(), 
					initial_lambdas.end(), 0.0);
  for(size_t i = 0; i < initial_lambdas.size(); i++)
    initial_lambdas[i] = initial_lambdas[i]/lambdas_sum;

  sort(initial_xs.begin(), initial_xs.end());
}

double
harris_by_newton(const bool VERBOSE,
		 const vector<double> &counts_hist,
		 const double tolerance,
		 const size_t max_iter,
		 const size_t depth){
  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];

  // depth corresponds to k in Harris
  vector<double> measure_moments;
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < depth; i++)
    measure_moments.push_back(exp(gsl_sf_lngamma(i + 3)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));


  // need to use an even number of coefficients
  // odd # gives upper bound
  if(depth % 2 == 1)
    measure_moments.pop_back();

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
      cerr << values_sum << endl;
      cerr << "lambdas = ";
      for(size_t i = 0; i < root_lambdas.size(); i++)
	cerr << root_lambdas[i] << ", ";
      cerr << endl;
    }
    for(size_t i = 0; i < root_xs.size(); i++)
      lower_bound += exp(log(counts_hist[1])
			 + log(root_lambdas[i])
			 - log(root_xs[i]));

    lower_bound += exp(log(counts_hist[1]) 
		       + log(root_lambdas.back()) 
		       -log(values_sum));
    lower_bound += 
      accumulate(counts_hist.begin(), counts_hist.end(), 0.0);
  }
  else
    cerr << "did not converge" << endl;
  

  return lower_bound;
}


///////////////////////////////////////////////////////
// Quadrature Methods to Estimate Library Size


double
golub_welsh_libsize(const bool VERBOSE,
		    const std::vector<double> &counts_hist,
		    const double tol,
		    const size_t max_iter,
		    const size_t n_points){

  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];

  size_t counts_before_first_zero = 1;
  while ((counts_before_first_zero < counts_hist.size()) &&
	  (counts_hist[counts_before_first_zero] > 0))
    ++counts_before_first_zero;
  if(2*n_points > counts_before_first_zero - 2)
    throw SMITHLABException("too many points for quadrature");
  
   
  // initialize moments, 0th moment is 1
  vector<double> measure_moments(1, 1.0);
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < 2*n_points; i++)
    measure_moments.push_back(exp(gsl_sf_lngamma(i + 3)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));

  vector<double> points, weights;
  golub_welsh_quadrature(VERBOSE, measure_moments, n_points,
			 tol,  max_iter, points, weights);

  // En_1 * int_0^\infty 1/x d \nu (x)
  double estimated_integral = 0.0;
  for(size_t i = 0; i < points.size(); i++)
    estimated_integral += counts_hist[1]*weights[i]/points[i];

  if(VERBOSE){
    cerr << "points = ";
    for(size_t i = 0; i < points.size(); i++)
      cerr << points[i] << ", ";
    cerr << endl;

    cerr << "weights = ";
    for(size_t i = 0; i < weights.size(); i++)
      cerr << weights[i] << ", ";
    cerr << endl;

    //  cerr << "estimated lib size = " << estimated_integral << endl;
  }

  return estimated_integral;
}

double
laguerre_modified_libsize(const bool VERBOSE,
			  const std::vector<double> &counts_hist,
			  const double tol,
			  const size_t max_iter,
			  const size_t n_points){
  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];

  size_t counts_before_first_zero = 1;
  while ((counts_before_first_zero < counts_hist.size()) &&
	  (counts_hist[counts_before_first_zero] > 0))
    ++counts_before_first_zero;
  if(2*n_points > counts_before_first_zero - 2)
    throw SMITHLABException("too many points for quadrature");
  
   
  // initialize moments, 0th moment is 1
  vector<double> measure_moments(1, 1.0);
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < 2*n_points; i++)
    measure_moments.push_back(exp(gsl_sf_lngamma(i + 3)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));

  vector<double> points, weights;
  laguerre_modified_quadrature(VERBOSE, measure_moments, n_points,
			       tol,  max_iter, points, weights);

  // En_1 * int_0^\infty 1/x d \nu (x)
  double estimated_integral = 0.0;
  for(size_t i = 0; i < points.size(); i++)
    estimated_integral += counts_hist[1]*weights[i]/points[i];

  if(VERBOSE){
    cerr << "points = ";
    for(size_t i = 0; i < points.size(); i++)
      cerr << points[i] << ", ";
    cerr << endl;

    cerr << "weights = ";
    for(size_t i = 0; i < weights.size(); i++)
      cerr << weights[i] << ", ";
    cerr << endl;

    //    cerr << "estimated lib size = " << estimated_integral << endl;
  }

  return estimated_integral;
}


double
chebyshev_libsize(const bool VERBOSE,
		  const std::vector<double> &counts_hist,
		  const double tol,
		  const size_t max_iter,
		  const size_t n_points){
  double values_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    values_sum += i*counts_hist[i];

  size_t counts_before_first_zero = 1;
  while ((counts_before_first_zero < counts_hist.size()) &&
	  (counts_hist[counts_before_first_zero] > 0))
    ++counts_before_first_zero;
  if(2*n_points > counts_before_first_zero - 2)
    throw SMITHLABException("too many points for quadrature");
  
   
  // initialize moments, 0th moment is 1
  vector<double> measure_moments(1, 1.0);
  // mu_r = (r + 1)! n_{r+1} / n_1
  for(size_t i = 0; i < 2*n_points; i++)
    measure_moments.push_back(exp(gsl_sf_lngamma(i + 3)
				  + log(counts_hist[i + 2])
				  - log(counts_hist[1])));

  vector<double> points, weights;
  chebyshev_quadrature(VERBOSE, measure_moments, n_points,
		       tol,  max_iter, points, weights);

  // En_1 * int_0^\infty 1/x d \nu (x)
  double estimated_integral = 0.0;
  for(size_t i = 0; i < points.size(); i++)
    estimated_integral += counts_hist[1]*weights[i]/points[i];

  if(VERBOSE){
    cerr << "points = ";
    for(size_t i = 0; i < points.size(); i++)
      cerr << points[i] << ", ";
    cerr << endl;

    cerr << "weights = ";
    for(size_t i = 0; i < weights.size(); i++)
      cerr << weights[i] << ", ";
    cerr << endl;

    //    cerr << "estimated lib size = " << estimated_integral << endl;
  }

  return estimated_integral;
}
