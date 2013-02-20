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

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_poly.h>



#include "smithlab_utils.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::cout;


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
poly_solve_gauss_quad(const std::vector<double> &moments,
		      const size_t n_points,
		      std::vector<double> &weights,
		      std::vector<double> &points){

  vector<double> moment_estimates(moments);
  moment_estimates.resize(2*n_points);

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
// Golub & Welsh quadrature

// Following Golub & Welsch 1968 (Sec 4)
static void
three_term_relation(const vector<double> &moments,
		    const size_t n_points,
		    vector<double> &a,
		    vector<double> &b){

  vector<double> moment_estimates(moments);
  moment_estimates.resize(2*n_points);

  // M is Hankel matrix of moments
  gsl_matrix *M = gsl_matrix_alloc(2*n_points, 2*n_points);
  for(size_t i = 0; i < n_points; i++){
    for(size_t j = 0; j < n_points; j++){
      gsl_matrix_set(M, i, j, moment_estimates[i + j]);
    }
  }

  // cholesky decomp on M
  // if M is not positive definite, error code GSL_EDOM should occur
  gsl_linalg_cholesky_decomp(M);

  // equations 4.3
  a.clear();
  for(size_t i = 0; i < n_points - 1; i++){
    double recurrence_val = gsl_matrix_get(M, i, i + 1)/gsl_matrix_get(M, i, i);
    if(i > 0)
      recurrence_val -= 
	gsl_matrix_get(M, i - 1, i)/gsl_matrix_get(M, i - 1, i - 1);
    a.push_back(recurrence_val);
  }

  b.clear();
  for(size_t i = 0; i < n_points - 2; i++)
    b.push_back(gsl_matrix_get(M, i + 1, i + 1)/gsl_matrix_get(M, i, i));

  gsl_matrix_free(M);  
}

// one iteration of QR: 
// following eq's 3.3 of Golub & Welsh
// make sure b is padded with an extra zero so b.size() = a.size()
static double
QR_iteration(vector<double> &a,
	     vector<double> &b,
	     vector<double> &z,
	     double lambda_guess){
  // error = b*b
  double return_error = 0.0;

  // initialize variables
  double b_bar = a[0] - lambda_guess;
  double d = b[0];
  double b_tilde = b[0];
  double z_bar = z[0];
  double sin_theta, cos_theta;
  vector<double> a_bar(a.size(), 0.0);
  a_bar[0] = a[0];

  for(size_t i = 0; i < a.size() - 1; i++){
    const double trig_denom = sqrt(d*d + b_bar*b_bar);
    sin_theta = d/trig_denom;
    cos_theta = b_bar/trig_denom;
    a[i] = a_bar*cos_theta*cos_theta + 2*b_tilde*cos_theta*sin_theta 
      + a[i + 1]*sin_theta*sin_theta;

    a_bar[i + 1] = a_bar[i]*sin_theta*sin_theta 
      - 2*b_tilde*cos_theta*sin_theta + a[i + 1]*cos_theta*cos_theta;
    b_bar = (a_bar[i] - a[i + 1])*sin_theta*cos_theta 
      + b_tilde*(sin_theta*sin_theta - cos_theta*cos_theta);
    z[i] = z_bar*cos_theta + z[i + 1]*sin_theta;

    if(i != 0){
      b[i - 1] = trig_denom;
      return_error + = b[i - 1]*b[i - 1];
    }

    if(i < a.size() - 2){
      b_tilde = -b[i + 1]*cos_theta;
      d = b[i + 1]*sin_theta; 
    }

    z_bar = z_bar*sin_theta - z[i + 1]*cos_theta;
  }

  //Last iteration? (not explicit in Golub & Welsh)
  const double trig_denom = sqrt(d*d + b_bar*b_bar);
  sin_theta = d/trig_denom;
  cos_theta = b_bar/trig_denom;
  a.back() = a_bar*cos_theta*cos_theta + 2*b_tilde*cos_theta*sin_theta;
  b.back() = trig_denom;
  z.back() = z_bar*cos_theta;
  return_error += b.back()*b.back();

  return return_error;
}


void
golub_welsh_quadrature(const vector<double> &moments,
		       const size_t n_points,
		       const double tol, const size_t max_iter,
		       vector<double> &points,
		       vector<double> &weights){
  // compute the 3-term recursion that generates the 
  // orthogonal polynomials
  vector<double> alpha, beta;
  three_term_relation(moments, n_points, alpha, beta);

  // add a zero to the back of beta so it's same size as alpha
  beta.push_back(0.0);

  vector<double> eigenvec(alpha.size(), 0.0);
  eigenvec[0] = 1.0;

  // can change lambda guess to affect convergence
  double lambda_guess = 0.0;

  // in QR, off-diagonals go to zero
  // use off diags for convergence
  double error = 0.0;
  for(size_t i = 0; i < beta.size(); i++)
    error += beta[i]*beta[i];
  size_t iter = 0;
  while(error >= tol && iter < max_iter){
    error = QR_iteration(alpha, beta, eigenvec, lambda_guess);
    iter++;
  }
  // eigenvalues are on diagonal of J, i.e. alpha
  points.swap(alpha);

  weights.swap(eigenvec);
  for(size_t i = 0; i < weights.size(); i++)
    weights[i] = weights[i]*weights[i];
}



////////////////////////////////////////////////////////////
//modified Chebyshev quadrature

// m[j] = 'th modified moment, v[j]=j'th moment
// generalized laguerre polynomial w/ a=1 : l_{j}(x)
// orthogonal to x e^{-x}
// l_j(x) = \sum_l=0^j 1/l! binom{1+j}{j-l} (-1)^l x^l
// m[j] = \sum_{l=0}^j 1/l! binom{1+j}{l-l} (-1)^l v[l]
static void
laguerre_modified_moments(const vector<double> &moments,
			  const size_t n_points,
			  vector<double> &modified_moments){
  modified_moments.resize(2*n_points, 0.0);
  for(size_t j = 0; j < modified_moments.size(); j++){
    for(size_t l = 0; l <= j; l++){
      modified_moments[j] += 
	pow(-1, l)*exp(- gsl_sf_lnfact(l) 
		       + gsl_sf_lnfact(j + 1)
		       - gsl_sf_lnfact(l + 1) 
		       - gsl_sf_lnfact(j - l))*moments[l]; 
    } 
  } 
}

// 3-term relation calculation by 
// modified Chabyshev algorithm
// Golub & Meurant (2010) pg 60 (bottom)
// a & b are the known 3-term relation
// of the modifying polynomials l_{j} (x),
// i.e. b_{j} l_{j+1} (x) = (x - a_j) l_j(x) - c_j l_{j-1} (x)
static void
modified3term_relation(const vector<double> &modified_moments,
		       const vector<double> &a,
		       const vector<double> &b,
		       const vector<double> &c,
		       const size_t n_points,
		       vector<double> &alpha,
		       vector<double> &nu){
  alpha.resize(n_points, 0.0);
  nu.resize(n_points, 0.0);

  vector< vector<double> > sigma(2*n_points, vector<double>(2*n_points, 0.0));
  // initialization
  alpha[0] = a[0] + modified_moments[0]/modified_moments[1];
  // sigma[-1][l] = 0
  for(size_t l = 0; l < 2*n_points, l++)
    sigma[0][l] = modified_moments[l];

  for(size_t k = 1; k < n_points; k++){
    for(size_t l = k; l < 2*n_points - k; l++){
      sigma[k][l] = b[l]*sigma[k-1][l+1] + (a[l] - alpha[k-1])*sigma[k-1][l]
	+ c[l]*sigma[k-1][l-1];
      if(k > 1)
	sigma[k][l] -= nu[k-2]*sigma[k-2][l];
    }
    alpha[k] = a[k] + b[k]*sigma[k][k+1]/sigma[k][k] - b[k-1]*sigma[k-1][k]/sigma[k-1][k-1];
    nu[k-1] = b[k-1]*sigma[k][k]/sigma[k-1][k-1];
  }  
}


void
laguerre_modified_quadrature(const vector<double> &moments,
			     const size_t n_points,
			     const double tol, const size_t max_iter,
			     vector<double> &points,
			     vector<double> &weights){
  // change of basis to laguerre polynomials
  vector<double> modified_moments;
  laguerre_modified_moments(moments, n_points, modified_moments);
 
  // b_{i+1} p_{i+1}(x) = (x - a_{i+1}) p_{i}(x) - c_i p_{i-2}(x)
  // -i l_i (x) = (x - 2*i) l_{i-1} (x) + i l_{i-2} (x)
  vector<double> a, b, c;
  for(size_t i = 1; i <= 2*n_points; i++){
    b.push_back(-i);
    a.push_back(2*i);
    c.push_back(-i); 
  }

  vector<double> alpha, beta;
  modified3term_relation(modified_moments, a, b, c,
			 n_points, alpha, beta);

  // add a zero to the back of beta so it's same size as alpha
  beta.push_back(0.0);

  vector<double> eigenvec(alpha.size(), 0.0);
  eigenvec[0] = 1.0;

  // can change lambda guess to affect convergence
  double lambda_guess = 0.0;

  // in QR, off-diagonals go to zero
  // use off diags for convergence
  double error = 0.0;
  for(size_t i = 0; i < beta.size(); i++)
    error += beta[i]*beta[i];
  size_t iter = 0;
  while(error >= tol && iter < max_iter){
    error = QR_iteration(alpha, beta, eigenvec, lambda_guess);
    iter++;
  }
  // eigenvalues are on diagonal of J, i.e. alpha
  points.swap(alpha);

  weights.swap(eigenvec);
  for(size_t i = 0; i < weights.size(); i++)
    weights[i] = weights[i]*weights[i];
}
