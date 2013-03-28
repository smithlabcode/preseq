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

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_poly.h>

#include "moment_sequence.hpp"

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cassert>

using std::string;
using std::vector;
using std::endl;
using std::max;
using std::cerr;



/////////////////////////////////////////////////////
// 3 term relations

// check 3 term recurrence to avoid non-positive elements
// truncate if non-positive element found
static void
check_three_term_relation(vector<double> &a,
			  vector<double> &b){

  // first entry is zero! Abort
  if(a[0] <= 0.0){
    a.clear();
    b.clear();
  }

  for(size_t i = 0; i < b.size(); i++){
    if(b[i] <= 0.0 || !finite(b[i])
       || a[i + 1] <= 0.0 || !finite(a[i + 1])){
      b.resize(i);
      a.resize(i + 1);
      break;
    }
  }
}

// check moment_sequence to avoid non-positive elements
// truncate if non-positive element found
static void
check_moment_sequence(vector<double> &obs_moms){
  if(obs_moms[0] <= 0.0 || !finite(obs_moms[0]))
     obs_moms.clear();

  for(size_t i = 1; i < obs_moms.size(); i++){
    if(obs_moms[i] <= 0.0 || !finite(obs_moms[i])){
      obs_moms.resize(i + 1);
      break;
    }
  }
}



/////////////////////////////////////////////////////
// Golub & Welsh quadrature

// Following Golub & Welsch 1968 (Sec 4)
void
MomentSequence::gw_three_term_calc(const bool VERBOSE,
				   const size_t n_points){

  const size_t num_points = 
    std::min(n_points, static_cast<size_t>(floor(moments.size()/2)) - 1);
  vector<double> moment_estimates(moments);
  moment_estimates.resize(2*num_points + 2);

  // M is Hankel matrix of moments
  gsl_matrix *M = gsl_matrix_alloc(num_points + 1, num_points + 1);
  for(size_t i = 0; i < num_points + 1; i++){
    for(size_t j = 0; j < num_points + 1; j++)
      gsl_matrix_set(M, i, j, moments[i + j]);
    
  }

  if(VERBOSE){
    cerr << "Hankel matrix = " << endl;
    for(size_t i = 0; i < num_points + 1; i++){
      for(size_t j  = 0; j < num_points + 1; j++)
	cerr << gsl_matrix_get(M, i, j) << ", ";
      cerr << endl;
    }
  }


  // cholesky decomp on M
  // if M is not positive definite, error code GSL_EDOM should occur
  gsl_linalg_cholesky_decomp(M);

  // equations 4.3
  vector<double> a, b;
  for(size_t i = 0; i < num_points; i++){
    double recurrence_val = gsl_matrix_get(M, i, i + 1)/gsl_matrix_get(M, i, i);
    if(i > 0)
      recurrence_val -= 
	gsl_matrix_get(M, i - 1, i)/gsl_matrix_get(M, i - 1, i - 1);
    a.push_back(recurrence_val);
  }

  for(size_t i = 0; i < num_points - 1; i++)
    b.push_back(gsl_matrix_get(M, i + 1, i + 1)/gsl_matrix_get(M, i, i));
  
  gsl_matrix_free(M);  

  // check recurrence
  check_three_term_relation(a, b);

  alpha = a;
  beta = b;
}



////////////////////////////////////////////////////////////
//modified Chebyshev quadrature


// 3-term relation calculation by 
// modified Chabyshev algorithm
// Golub & Meurant (2010) pg 60 (bottom)
// a & b are the known 3-term relation
// of the modifying polynomials l_{j} (x),
// i.e. b_{j} l_{j+1} (x) = (x - a_j) l_j(x) - c_j l_{j-1} (x)
void
MomentSequence::modified_Chebyshev(const bool VERBOSE,
				   const size_t n_points,
				   const vector<double> &mod_alpha,
				   const vector<double> &mod_beta,
				   const vector<double> &modified_moments){
  vector<double> a(n_points, 0.0);
  vector<double> b(n_points - 1, 0.0);

  vector< vector<double> > sigma(2*n_points, vector<double>(2*n_points, 0.0));
  // initialization
  a[0] = mod_alpha[0] + modified_moments[1]/modified_moments[0];
  // sigma[-1][l] = 0
  for(size_t l = 0; l < 2*n_points; l++)
    sigma[0][l] = modified_moments[l];

  for(size_t k = 1; k <= n_points; k++){
    for(size_t l = k; l < 2*n_points - k; l++){
      sigma[k][l] = sigma[k-1][l+1] + (mod_alpha[l] - a[k-1])*sigma[k-1][l]
	+ mod_beta[l-1]*sigma[k-1][l-1];
      if(k > 1)
	sigma[k][l] -= b[k-2]*sigma[k-2][l];
    }
    if(k != n_points){
      a[k] = mod_alpha[k] + sigma[k][k+1]/sigma[k][k] - sigma[k-1][k]/sigma[k-1][k-1];
      b[k-1] = sigma[k][k]/sigma[k-1][k-1];
    }
  }  

  check_three_term_relation(a, b);

  // See Gautschi pgs 10-13,
  // the nu here is the square of the off-diagonal
  // of the Jacobi matrix
  for(size_t i = 0; i < b.size(); i++)
    b[i] = sqrt(b[i]);

  alpha = a;
  beta = b;
}

void
MomentSequence::unmodified_Chebyshev(const bool VERBOSE){

  const size_t n_points = static_cast<size_t>(floor(moments.size()/2));
  vector<double> a(n_points, 0.0);
  vector<double> b(n_points - 1, 0.0);

  vector< vector<double> > sigma(2*n_points, vector<double>(2*n_points, 0.0));
  // initialization
  a[0] = moments[1]/moments[0];
  // sigma[-1][l] = 0
  for(size_t l = 0; l < 2*n_points; l++)
    sigma[0][l] = moments[l];

  for(size_t k = 1; k <= n_points; k++){
    for(size_t l = k; l < 2*n_points - k; l++){
      sigma[k][l] = sigma[k-1][l+1] - a[k-1]*sigma[k-1][l];
      if(k > 1)
	sigma[k][l] -= b[k-2]*sigma[k-2][l];
    }
    if(k != n_points){
      a[k] = sigma[k][k+1]/sigma[k][k] - sigma[k-1][k]/sigma[k-1][k-1];
      b[k-1] = sigma[k][k]/sigma[k-1][k-1];
    }
  }  

  if(VERBOSE){
    cerr << "3-term relations:" << endl;
    cerr << "alpha = ";
    for(size_t i = 0; i < a.size(); i++)
      cerr << a[i] << ", ";
    cerr << endl;
    cerr << "beta = ";
    for(size_t i = 0; i < b.size(); i++)
      cerr << b[i] << ", ";
    cerr << endl;
  }

  check_three_term_relation(a, b);

  // See Gautschi pgs 10-13,
  // the nu here is the square of the off-diagonal
  // of the Jacobi matrix
  for(size_t i = 0; i < b.size(); i++)
    b[i] = sqrt(b[i]);

  alpha = a;
  beta = b;
}


////////////////////////////////////////////////////
// Constructor

MomentSequence::MomentSequence(const vector<double> &obs_moms) :
  moments(obs_moms) {
  vector<double> holding_moms(moments);
  // make sure the moments are all positive
  check_moment_sequence(holding_moms);
  moments = holding_moms;

  // calculate 3-term recurrence
  unmodified_Chebyshev(true);
}


/////////////////////////////////////////////////////
// Quadrature rule calculations

static void
solve_linear_system(const vector<vector<double> > &A_in,
		    const vector<double> &b_in,
		    vector<double> &x_out) {
  const size_t rows = A_in.size();
  const size_t cols = A_in.front().size();
  
  assert(b_in.size() >= rows);
  
  gsl_matrix *A = gsl_matrix_alloc(rows, cols);
  
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j)
      gsl_matrix_set(A, i, j, A_in[i][j]);
  
  gsl_vector *b = gsl_vector_alloc(rows);
  
  for (size_t i = 0; i < rows; ++i)
    gsl_vector_set(b, i, -1*b_in[rows + i]);
  
  gsl_vector *x = gsl_vector_alloc(cols);
  gsl_linalg_HH_solve(A, b, x);
  //gsl_linalg_LU_solve(A, b, x);

  for (size_t i = 0; i < cols; ++i)
    x_out.push_back(gsl_vector_get(x, i));
  
  gsl_matrix_free(A);
  gsl_vector_free(b);
  gsl_vector_free(x);
}   

void
MomentSequence::poly_solve_gauss_quad(const size_t n_points,
				      vector<double> &weights,
				      vector<double> &points){

  // make sure n_points can be used
  vector<double> moment_estimates(moments);
  moment_estimates.resize((n_points < moments.size()/2) ? 2*n_points : moments.size());

  vector<vector<double> > matrix(n_points, vector<double>(n_points, 0.0));
  for (size_t i = 0; i < n_points; ++i) {
    for (size_t j = 0; j < n_points; ++j) {
      matrix[i][j] = moment_estimates[j + i];
    }
  }

  vector<double> c;
  solve_linear_system(matrix, moment_estimates, c);
  c.push_back(1);

   gsl_poly_complex_workspace *w =
      gsl_poly_complex_workspace_alloc(n_points + 1);
    
   vector<double> x(2*(n_points + 1), 0.0);
   gsl_poly_complex_solve(&c[0], c.size(), w, &x[0]);

   weights.swap(c);
   points.swap(x);
   gsl_poly_complex_workspace_free (w);
}


/////////////////////////////////////////////////////
// Quadrature Methods

// one iteration of QR: 
// following eq's 3.3 of Golub & Welsh
// one iteration is Z_N-1*Z_N-2*...*Z_1*X*Z_1*...*Z_N-1
// Z_j is givens matrix to zero out the j+1,j'th element of X
static void
QRiteration(vector<double> &alpha,
	    vector<double> &beta,
	    vector<double> &weights){
  // initialize variables
  vector<double> sin_theta(alpha.size(), 0.0);
  vector<double> cos_theta(alpha.size(), 0.0);

  vector<double> a(alpha.size(), 0.0);
  vector<double> a_bar(alpha.size(), 0.0);
  a_bar[0] = alpha[0];

  vector<double> b(beta);
  vector<double> b_bar(alpha.size(), 0.0);
  b_bar[0] = alpha[0];
  vector<double> b_tilde(alpha.size(), 0.0);
  b_tilde[0] = beta[0];

  vector<double> d(alpha.size(), 0.0);
  d[0] = beta[0];

  vector<double> z(weights);
  vector<double> z_bar(weights.size(), 0.0);
  z_bar[0] = z[0];


  for(size_t j = 0; j < alpha.size() - 1; j++){
    // for d and b_bar, j here is j-1 in G&W 
    sin_theta[j] = d[j]/sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);
    cos_theta[j] = b_bar[j]/sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);

    a[j] = a_bar[j]*cos_theta[j]*cos_theta[j]
      + 2*b_tilde[j]*cos_theta[j]*sin_theta[j] 
      + alpha[j+1]*sin_theta[j]*sin_theta[j];

    a_bar[j+1] = a_bar[j]*sin_theta[j]*sin_theta[j]
      - 2*b_tilde[j]*cos_theta[j]*sin_theta[j] 
      + alpha[j+1]*cos_theta[j]*cos_theta[j];

    if(j != 0)
      b[j-1] = sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);

    b_bar[j+1] = (a_bar[j] - alpha[j+1])*sin_theta[j]*cos_theta[j]
      + b_tilde[j]*(sin_theta[j]*sin_theta[j] - cos_theta[j]*cos_theta[j]);

    b_tilde[j+1] = -beta[j+1]*cos_theta[j];
    
    d[j+1] = beta[j+1]*sin_theta[j];
 
    z[j] = z_bar[j]*cos_theta[j] + weights[j+1]*sin_theta[j];

    z_bar[j+1] = z_bar[j]*sin_theta[j] - weights[j+1]*cos_theta[j];
  }

// last entries set equal to final "holding" values
  a.back() = a_bar.back();
  b.back() = b_bar.back();
  z.back() = z_bar.back();

  alpha.swap(a);
  beta.swap(b);
  weights.swap(z);
}

static bool
check_positivity(const vector<double> &points){
  for(size_t i = 0; i < points.size(); i++)
    if(points[i] <= 0.0 || !finite(points[i]))
      return false;

  return true;
}


void
MomentSequence::QR_quadrature_rules(const bool VERBOSE,
				    const size_t n_points,
				    const double tol, 
				    const size_t max_iter,
				    vector<double> &points,
				    vector<double> &weights){

  // make sure that points.size() will be less than n_points
  vector<double> a(alpha);
  a.resize((n_points < alpha.size()) ? n_points : alpha.size());
  vector<double> b(beta);
  b.resize((n_points - 1 < beta.size()) ? n_points - 1 : beta.size());

  bool POSITIVE_POINTS = false;

  while(!(POSITIVE_POINTS) && a.size() > 0){
    vector<double> eigenvec(a.size(), 0.0);
    eigenvec[0] = 1.0;
    vector<double> eigenvals(a);
    vector<double> qr_beta(b);
  // in QR, off-diagonals go to zero
  // use off diags for convergence
    double error = 0.0;
    for(size_t i = 0; i < qr_beta.size(); i++)
      error += qr_beta[i]*qr_beta[i];
    size_t iter = 0;
    while(iter < max_iter){
      QRiteration(eigenvals, qr_beta, eigenvec);

      error = 0.0;
      for(size_t i = 0; i < qr_beta.size(); i++)
	error += qr_beta[i]*qr_beta[i];
      iter++;
    }
  // eigenvalues are on diagonal of J
    POSITIVE_POINTS = check_positivity(eigenvals);

    if(VERBOSE){
      cerr << "POINTS = ";
      for(size_t i = 0; i < eigenvals.size(); i++)
	cerr << eigenvals[i] << ", ";
      cerr << endl;
    }

    if(POSITIVE_POINTS){
      points.swap(eigenvals);
      weights.swap(eigenvec);
    }
    else{
      a.pop_back();
      b.pop_back();
    }

  }

  for(size_t i = 0; i < weights.size(); i++)
    weights[i] = weights[i]*weights[i];
}



