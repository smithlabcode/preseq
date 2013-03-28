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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

#include "smithlab_utils.hpp"

using std::string;
using std::vector;
using std::max;
using std::cerr;
using std::endl;
using std::numeric_limits;


using smithlab::log_sum_log_vec;


// to pass parameters through
struct parameters{
  vector<double> lambdas;
  vector<double> xs;
  vector<double> moments;
  double vals_sum;
};

// lambdas_xs is lambdas followed by x's
// last x = N
// eqns: sum_i lambda_i = 1
// sum_i lambda_i x_i = mu_1
// sum_i lambda_i x_i^2 = mu_2
// ....
// sum_i lambda_i x_i^k = mu_k
static void
system_eqns(const parameters &params,
	    vector<double> &output){
  vector<double> moms, lambs, x;
  moms = params.moments;
  lambs = params.lambdas;
  x = params.xs;
  output.resize(moms.size(), 0.0);

  for(size_t i = 0; i < moms.size(); i++){
    double eqn_value = 0.0;
    for(size_t j = 0; j < x.size(); j++)
      eqn_value += lambs[j]*pow(x[j], i);
    output[i] = eqn_value - moms[i];
  }
}

static void 
jacobian(const parameters &params,
	 vector< vector<double> > &jacob){
  vector<double> lambs, x;
  lambs = params.lambdas;
  x = params.xs;

  const size_t variable_dim = lambs.size() + x.size();
  vector< vector<double> > computed_jacobian(variable_dim,
					     vector<double>(variable_dim, 0.0));

  // derivatives of sum_i lambda_i = 1
  for(size_t j = 0; j < lambs.size(); j++)
    computed_jacobian[0][j] = 1.0;

  // dervatives of sum_i lambda_i x_i^k = mu_k for lambdas
  for(size_t k = 1; k < computed_jacobian.size(); k++)
    for(size_t j = 0; j < lambs.size(); j++)
      computed_jacobian[k][j] = pow(x[j], k);

  // derivative of sum_i lambda_i x_i^k = mu_k for xs
  for(size_t k = 1; k < computed_jacobian.size(); k++)
    for(size_t i = 0; i < x.size(); i++)
      computed_jacobian[k][lambs.size() + i] = 
	k*lambs[i]*pow(x[i], k - 1);

  jacob.swap(computed_jacobian);
}

// converged if sum residuals^2 < tol
static inline bool
check_convergence(const vector<double> &current_values,
		  const double tolerance){

  vector<double> log_residuals;
  for(size_t i = 0; i < current_values.size(); i++)
    log_residuals.push_back(2*log(fabs(current_values[i])));

  if(exp(log_sum_log_vec(log_residuals, 
			 log_residuals.size())) < tolerance)
    return true;

  return false;
}

// check the x's to make sure all are strictly positive,
// as guaranteed by Harris
static inline bool
check_positive(const vector<double> &x){
  bool POS_X = true;
  for(size_t i = 0; i < x.size(); i++)
    POS_X = POS_X && (x[i] > 0.0);

  return POS_X;
}


static inline bool
check_finite(const vector<double> &x){
  bool ALL_FINITE = true;
  for(size_t i = 0; i < x.size(); i++)
    ALL_FINITE = ALL_FINITE && finite(x[i]);

  return ALL_FINITE;
}

static inline bool
check_xs(const vector<double> &xs,
	 const double N){
  // Check if they are increasing and less than N
  if(xs[0] > N)
    return false;
  if(xs.size() > 1){
    for(size_t i = 1; i < xs.size(); i++)
      if(xs[i] < xs[i - 1] || xs[i] > N)
	return false;
  }

  return true;
}

static inline bool
smaller_abs(const double a, const double b){
  return fabs(a) < fabs(b);
}

// we want the iteration to move in the same direction,
// but we want no variable to be negative.
static inline void
calculate_nonneg_delta(const parameters &params,
		       const vector<double> &current_lambdas,
		       const vector<double> &current_xs,
		       vector<double> &proposed_delta){

  const double upper_bound = params.vals_sum;

  // use L_infinity norm
  vector<double> full_vec(current_lambdas);
  full_vec.insert(full_vec.end(), current_xs.begin(), current_xs.end());
  assert(full_vec.size() == proposed_delta.size());

  /*
  cerr << "full_vec = ";
  for(size_t i = 0; i < full_vec.size(); i++)
    cerr << full_vec[i] << ", ";
  cerr << endl;
  */

  // multiplier = smallest y s.t. x_i + y*delta_x_i = 0
  double lower_multiplier = 
    ( (full_vec[0] + proposed_delta[0] <= 0.0) ? 
      -full_vec[0]/proposed_delta[0] : numeric_limits<double>::max() );
  for(size_t i = 1; i < full_vec.size(); i++){
    if(full_vec[i] + proposed_delta[i] <= 0.0){
      //    cerr << "indx " << i << " has neg new pos" << endl;
      lower_multiplier = std::min(lower_multiplier, -full_vec[i]/proposed_delta[i]);
    }
  }
  //go only half the distance to origin
  lower_multiplier = lower_multiplier/2.0;

  double upper_multiplier = 
    ( (full_vec[0] + proposed_delta[0] >= upper_bound) ? 
      (upper_bound - full_vec[0])/proposed_delta[0] : numeric_limits<double>::max() );

  for(size_t i = 1; i < full_vec.size(); i++){
    if(full_vec[i] + proposed_delta[i] >= upper_bound){
      //    cerr << "indx " << i << " has neg new pos" << endl;
      upper_multiplier = std::min(upper_multiplier, 
				  (upper_bound - full_vec[i])/proposed_delta[i]);
    }
  }
  //go only half the distance to upper boundary
  upper_multiplier = upper_multiplier/2.0;

  // some new val is negative, 
  // set new delta in same direction as proposed,
  // make sure all new variables are positive
  // set new delta so that the largest neg new val
  // is halved
  if(lower_multiplier < 1.0 || upper_multiplier < 1.0){
    double multiplier = std::min(lower_multiplier, upper_multiplier);
    for(size_t i = 0; i < proposed_delta.size(); i++)
      proposed_delta[i] = proposed_delta[i]*multiplier;

    /*
    cerr << "new_delta : ";
    for(size_t i = 0; i < proposed_delta.size(); i++)
      cerr << proposed_delta[i] << ", ";
    cerr << endl;
    */

  }
  // else the proposed_delta is acceptable, no new
  // vals are neg, do nothing
}

// find new guess on zero by householder trans of J*delta_x = - f
// return false if det(J) = 0
static bool
iterate_newton(const vector<double> &current_values,
	       const vector< vector<double> > &jacob,
	       const double tolerance,
	       parameters &params){
  vector<double> current_lambdas, current_xs;
  current_lambdas = params.lambdas;
  current_xs = params.xs;

  gsl_vector *neg_f_vals = gsl_vector_alloc(current_values.size());
  for(size_t i = 0; i < current_values.size(); i++)
    gsl_vector_set(neg_f_vals, i, -current_values[i]);

  gsl_matrix *LU = gsl_matrix_alloc(jacob.size(), jacob[0].size());
  gsl_matrix *QR = gsl_matrix_alloc(jacob.size(), jacob[0].size());
  for(size_t i = 0; i < jacob.size(); i++){
    for(size_t k = 0; k < jacob[i].size(); k++){
      gsl_matrix_set(QR, i, k, jacob[i][k]);
      gsl_matrix_set(LU, i, k, jacob[i][k]);
    }
  }

  // delta_x = (x_n+1 - x_n) i.e. change in x at current iteration
  gsl_vector *direction = gsl_vector_calloc(current_lambdas.size() +  
					    current_xs.size());

  gsl_permutation *P = gsl_permutation_alloc(jacob.size());
  int signum_P = 0;
  gsl_vector *tau = gsl_vector_calloc(jacob.size());

  gsl_linalg_LU_decomp(LU, P, &signum_P);
  gsl_linalg_QR_decomp(QR, tau);

  double log_det = gsl_linalg_LU_lndet(LU);

  if(finite(fabs(log_det))){

    gsl_linalg_LU_solve(LU, P, neg_f_vals, direction);
    //gsl_linalg_QR_solve(QR, tau, neg_f_vals, direction);

    vector<double> delta_x(current_lambdas.size() 
			   + current_xs.size(), 0.0);
    for(size_t i = 0; i < delta_x.size(); i++)
      delta_x[i] = gsl_vector_get(direction, i);

    calculate_nonneg_delta(params, current_lambdas, 
			   current_xs, delta_x);



    vector<double> new_lambdas = current_lambdas;
    for(size_t i = 0; i < current_lambdas.size(); i++)
      new_lambdas[i] += delta_x[i];

    // should we normalize lambdas?
    const double lambdas_sum = accumulate(new_lambdas.begin(), new_lambdas.end(), 0.0);
    for(size_t i = 0; i < new_lambdas.size(); i++)
      new_lambdas[i] = new_lambdas[i]/lambdas_sum;
    params.lambdas = new_lambdas;

    vector<double> new_xs = current_xs;
    for(size_t i = 0; i < current_xs.size(); i++)
      new_xs[i] += delta_x[current_lambdas.size() + i];
    
    // check new points, exit with failure if they don't check
    //  if(!check_xs(new_xs, params.vals_sum))
    //    return false;

    params.xs = new_xs;

    vector<double> new_values;
    system_eqns(params, new_values);
    /*
    cerr << "new lambdas: ";
    for(size_t i = 0; i < new_lambdas.size(); i++)
      cerr << new_lambdas[i] << ", ";
    cerr << endl;
    cerr << "new xs     : ";
    for(size_t i = 0; i < new_xs.size(); i++)
      cerr << new_xs[i] << ", ";
    cerr << endl;
    cerr << "new values : ";
    for(size_t i = 0 ; i < new_values.size(); i++)
      cerr << new_values[i] << ", ";
    cerr << endl;
    */
    

    gsl_matrix_free(LU);
    gsl_matrix_free(QR);
    gsl_permutation_free(P);
    gsl_vector_free(neg_f_vals);
    gsl_vector_free(direction);
    gsl_vector_free(tau);
    return true;

  }
  /*
  else{
    cerr << "determinant too small, log det = " << log_det << ", restarting from rand pos" << endl;
    cerr << "JACOBIAN = " << endl;
    for(size_t i = 0; i < jacob.size(); i++){
      for(size_t j = 0; j < jacob[i].size(); j++)
	cerr << jacob[i][j] << ", ";
      cerr << endl;
    }
    cerr << "CURRENT LAMBDAS = " << endl;
    for(size_t i = 0 ; i < params.lambdas.size(); i++)
      cerr << params.lambdas[i] << ", ";
    cerr << endl;
    cerr << "CURRENT Xs = " << endl;
    for(size_t i = 0; i < params.xs.size(); i++)
      cerr << params.xs[i] << ", ";
    cerr << endl;
    cerr << "CURRENT VALUES = " << endl;
    for(size_t i = 0; i < current_values.size(); i++)
      cerr << current_values[i] << ", ";
    cerr << endl;
    throw SMITHLABException("Fuck");
  }
  */
  gsl_matrix_free(LU);
  gsl_matrix_free(QR);
  gsl_permutation_free(P);
  gsl_vector_free(neg_f_vals);
  gsl_vector_free(direction);
  gsl_vector_free(tau);

  // if condition number too large, 
  // exit and flag iteration as unsuccesful
  return false;
}

static bool 
full_iteration_newton(const bool VERBOSE,
		      const double tolerance,
		      const size_t max_iter,
		      const size_t degrees_freedom,
		      parameters &params){
  vector<double> current_func_vals(degrees_freedom, 0.0);
  vector< vector<double> > J(degrees_freedom, vector<double>(degrees_freedom, 0.0));

  vector<double> old_lambdas, old_xs;
  system_eqns(params, current_func_vals);
  size_t indx = 0;
  bool CONVERGED = false;
  bool ITERATION_SUCCESS = true;
  bool CONTINUE = true;
  do{
    jacobian(params, J);

    ITERATION_SUCCESS = iterate_newton(current_func_vals, J, tolerance, params);

    system_eqns(params, current_func_vals);
    /*
    if(VERBOSE){
      cerr << "new lambdas : ";
      for(size_t i = 0; i < params.lambdas.size(); i++)
	cerr << params.lambdas[i] << ", ";
      cerr << endl;
      cerr << "new xs : ";
      for(size_t i = 0; i < params.xs.size(); i++)
	cerr << params.xs[i] << ", ";
      cerr << endl;
    }
    */

    CONVERGED = check_convergence(current_func_vals, tolerance);
    indx++;	
    if(!ITERATION_SUCCESS)
      cerr << "iteration unsuccessful" << endl;

  //exit loop if iteration unsuccessful(det = 0), 
  //indx reaches max_iter, or convergence 
    CONTINUE = ITERATION_SUCCESS && (indx < max_iter) && !CONVERGED;
  }while(CONTINUE);

  //  cerr << "new starting point" << endl;

  //
  if(CONVERGED && ITERATION_SUCCESS){
    //check all estimates are positive 
    return check_positive(params.xs) && check_positive(params.lambdas);
  }
  
  return false;
}


bool
modified_newtons_method(const bool VERBOSE,
			const vector<double> &initial_lambdas,
			const vector<double> &initial_xs,
			const vector<double> &in_moments,
			const double values_sum,
			const double tolerance, const size_t max_iter,
			vector<double> &root_lambdas,
			vector<double> &root_xs){

  const size_t degrees_freedom = initial_lambdas.size() + initial_xs.size();
  // set params
  parameters params;
  params.lambdas = initial_lambdas;
  params.xs = initial_xs;
  params.moments = in_moments;
  params.vals_sum = values_sum;

  bool FULL_ITER_SUCCESS = true;
  FULL_ITER_SUCCESS = full_iteration_newton(VERBOSE, tolerance, max_iter,
					    degrees_freedom, params);

  // make sure all values are finite
  FULL_ITER_SUCCESS =
    FULL_ITER_SUCCESS && check_finite(params.lambdas) && check_finite(params.xs);
  if(!check_finite(params.lambdas) || !check_finite(params.xs))
    cerr << "values are not finite, start over." << endl;

  if(FULL_ITER_SUCCESS){
    root_lambdas = params.lambdas;
    root_xs = params.xs;
    return true;
  }
  else{
    root_lambdas.clear();
    root_xs.clear();
  }

  // return false if iteration unsuccesful
  return false;
}
