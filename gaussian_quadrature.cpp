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
  cerr << "Hankel matrix:" << endl;
  gsl_matrix *M = gsl_matrix_alloc(n_points, n_points);
  for(size_t i = 0; i < n_points; i++){
    for(size_t j = 0; j < n_points; j++){
      gsl_matrix_set(M, i, j, moment_estimates[i + j]);
      cerr << gsl_matrix_get(M, i, j) << ", ";
    }
    cerr << endl;
  }


  // cholesky decomp on M
  // if M is not positive definite, error code GSL_EDOM should occur
  gsl_linalg_cholesky_decomp(M);

  cerr << "Cholesky decomp:" << endl;
  for(size_t i = 0; i < n_points; i++){
    for(size_t j = 0; j < n_points; j++){
      cerr << gsl_matrix_get(M, i, j) << ", ";
    }
    cerr << endl;
  }

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


void
golub_welsh_quadrature(const bool VERBOSE,
		       const vector<double> &moments,
		       const size_t n_points,
		       const double tol, const size_t max_iter,
		       vector<double> &points,
		       vector<double> &weights){
  // compute the 3-term recursion that generates the 
  // orthogonal polynomials
  vector<double> alpha, beta;
  three_term_relation(moments, n_points, alpha, beta);

  if(VERBOSE){
    cerr << "moments = ";
    for(size_t i = 0; i < 2*n_points; i++)
      cerr << moments[i] << ", ";
    cerr << endl;

    cerr << "alpha = ";
    for(size_t i = 0; i < alpha.size(); i++)
      cerr << alpha[i] << ", ";
    cerr << endl;

    cerr << "beta = ";
    for(size_t i = 0; i < beta.size(); i++)
      cerr << beta[i] << ", ";
    cerr << endl;
  }


  vector<double> eigenvec(alpha.size(), 0.0);
  eigenvec[0] = 1.0;

  // in QR, off-diagonals go to zero
  // use off diags for convergence
  double error = 0.0;
  for(size_t i = 0; i < beta.size(); i++)
    error += beta[i]*beta[i];
  size_t iter = 0;
  while(error >= tol && iter < max_iter){
    QRiteration(alpha, beta, eigenvec);

    error = 0.0;
    for(size_t i = 0; i < beta.size(); i++)
      error += beta[i]*beta[i];
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
// monic generalized laguerre polynomial w/ a=1 : l_{j}(x)
// orthogonal to x e^{-x}
// l_j(x) = \sum_l=0^j j!/l! binom{1+j}{j-l} (-1)^{l+j} x^l
// m[j] = \sum_{l=0}^j j!/l! binom{1+j}{j-l} (-1)^{l+j} v[l]
static void
laguerre_modified_moments(const vector<double> &moments,
			  const size_t n_points,
			  vector<double> &modified_moments){
  modified_moments.resize(2*n_points, 0.0);
  for(size_t j = 0; j < modified_moments.size(); j++){
    for(size_t l = 0; l <= j; l++){
      modified_moments[j] += 
	exp(gsl_sf_lnfact(j) - gsl_sf_lnfact(l)
	    + gsl_sf_lnfact(1 + j) - gsl_sf_lnfact(j - l)
	    - gsl_sf_lnfact(l + 1) + log(moments[l]))*pow(-1, l + j); 
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
		       const vector<double> &c,
		       const size_t n_points,
		       vector<double> &alpha,
		       vector<double> &nu){
  alpha.resize(n_points, 0.0);
  nu.resize(n_points - 1, 0.0);

  vector< vector<double> > sigma(2*n_points, vector<double>(2*n_points, 0.0));
  // initialization
  alpha[0] = a[0] + modified_moments[1]/modified_moments[0];
  // sigma[-1][l] = 0
  for(size_t l = 0; l < 2*n_points; l++)
    sigma[0][l] = modified_moments[l];

  for(size_t k = 1; k <= n_points; k++){
    for(size_t l = k; l < 2*n_points - k + 1; l++){
      sigma[k][l] = sigma[k-1][l+1] + (a[l] - alpha[k-1])*sigma[k-1][l]
	+ c[l]*sigma[k-1][l-1];
      if(k > 1)
	sigma[k][l] -= nu[k-2]*sigma[k-2][l];
    }
    alpha[k] = a[k] + sigma[k][k+1]/sigma[k][k] - sigma[k-1][k]/sigma[k-1][k-1];
    nu[k-1] = sigma[k][k]/sigma[k-1][k-1];
  }  

  // See Gautschi pgs 10-13,
  // the nu here is the square of the off-diagonal
  // of the Jacobi matrix
  for(size_t i = 0; i < nu.size(); i++)
    nu[i] = sqrt(nu[i]);
}


void
laguerre_modified_quadrature(const bool VERBOSE,
			     const vector<double> &moments,
			     const size_t n_points,
			     const double tol, const size_t max_iter,
			     vector<double> &points,
			     vector<double> &weights){
  // change of basis to laguerre polynomials
  vector<double> modified_moments;
  laguerre_modified_moments(moments, n_points-1, modified_moments);

  if(VERBOSE){
    cerr << "ORIGINAL MOMENTS = ";
    for(size_t i = 0; i < moments.size(); i++)
      cerr << moments[i] << ", ";
    cerr << endl;

    cerr << "MODIFIED MOMENTS = ";
    for(size_t i = 0; i < modified_moments.size(); i++)
      cerr << modified_moments[i] << ", ";
    cerr << endl;
  }
 
  //  p_{i+1}(x) = (x - a_{i+1}) p_{i}(x) - b_i p_{i-2}(x)
  // l_i (x) = (x - 2*(i+1)) l_{i-1} (x) - i*(i+1) l_{i-2} (x)
  vector<double> a, b;
  for(size_t i = 0; i < 2*(n_points-1); i++){
    b.push_back(i*(i + 1.0));
    a.push_back(2.0*(i+1));
  }

  vector<double> alpha, beta;
  modified3term_relation(modified_moments, a, b,
			 n_points-1, alpha, beta);

  if(VERBOSE){
    cerr << "ORIGINAL 3-TERM RELATION:" << endl;
    cerr << "a = ";
    for(size_t i = 0; i < a.size(); i++)
      cerr << a[i] << ", ";
    cerr << endl;
    cerr << "b = ";
    for(size_t i = 0; i < b.size(); i++)
      cerr << b[i] << ", ";
    cerr << endl;

    cerr << "ESTIMATED 3=TERM RELATION:" << endl;
    cerr << "alpha = ";
    for(size_t i = 0; i < alpha.size(); i++)
      cerr << alpha[i] << ", ";
    cerr << endl;
    cerr << "beta = ";
    for(size_t i = 0; i < beta.size(); i++)
      cerr << beta[i] << ", ";
    cerr << endl; 
  }
  vector<double> eigenvec(alpha.size(), 0.0);
  eigenvec[0] = 1.0;

  // in QR, off-diagonals go to zero
  // use off diags for convergence
  double error = 0.0;
  for(size_t i = 0; i < beta.size(); i++)
    error += beta[i]*beta[i];
  size_t iter = 0;
  while(error >= tol && iter < max_iter){
    QRiteration(alpha, beta, eigenvec);

    error = 0.0;
    for(size_t i = 0; i < beta.size(); i++)
      error += beta[i]*beta[i];
    iter++;
  }
  // eigenvalues are on diagonal of J, i.e. alpha
  points.swap(alpha);

  weights.swap(eigenvec);
  for(size_t i = 0; i < weights.size(); i++)
    weights[i] = weights[i]*weights[i];

}
