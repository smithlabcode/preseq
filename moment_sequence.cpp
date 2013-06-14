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
#include <gsl/gsl_randist.h>

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
using std::setprecision;



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

// un-normalized 3 term recurrence
void
MomentSequence::full_3term_recurrence(const bool VERBOSE,
				      vector<double> &full_alpha,
				      vector<double> &full_beta){

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

  full_alpha.swap(a);
  full_beta.swap(b);
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
  unmodified_Chebyshev(false);
}

//////////////////////////////////////////////////
// mean 3 term recurrence by bootstrapping histogram

static void
resample_hist(const gsl_rng *rng, const vector<double> &vals_hist,
	      const double total_sampled_reads,
	      double expected_sample_size,
	      vector<double> &sample_hist) {
  
  const size_t hist_size = vals_hist.size();
  const double vals_mean = total_sampled_reads/expected_sample_size;
  
  sample_hist = vector<double>(hist_size, 0.0);
  vector<unsigned int> curr_sample(hist_size);
  double remaining = total_sampled_reads;
  
  while (remaining > 0) {
    
    // get a new sample
    expected_sample_size = max(1.0, (remaining/vals_mean)/2.0);
    gsl_ran_multinomial(rng, hist_size, 
			static_cast<unsigned int>(expected_sample_size),
			&vals_hist.front(), &curr_sample.front());
    
    // see how much we got
    double inc = 0.0;
    for (size_t i = 0; i < hist_size; ++i)
      inc += i*curr_sample[i];
    
    // only add to histogram if sampled reads < remaining reads
    if (inc <= remaining) {
      for (size_t i = 0; i < hist_size; i++)
	sample_hist[i] += static_cast<double>(curr_sample[i]);
      // update the amount we still need to get
      remaining -= inc;
    }
  }
}

void
MomentSequence::set_mean_3term_recurrence(const bool VERBOSE, 
					  const vector<double> &counts_hist,
					  const size_t bootstraps,
					  const size_t n_points){
  double vals_sum = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    vals_sum += i*counts_hist[i];

  const double observed_distinct = accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, rand()); 

  vector<double> mean_alpha(n_points, 0.0);
  vector<double> mean_beta(n_points - 1, 0.0);

  for(size_t iter = 0; iter < bootstraps; iter++){

    vector<double> sample_hist;
    resample_hist(rng, counts_hist, vals_sum, observed_distinct, sample_hist);

    // construct sample moments
    vector<double> sample_moments;
  // mu_r = (r + 1)! n_{r+1} / n_1
    size_t indx = 1;
    while(sample_hist[indx] > 0  && indx <= std::min(2*n_points, sample_hist.size())){
      sample_moments.push_back(exp(gsl_sf_lnfact(indx)
				   + log(sample_hist[indx])
				   - log(sample_hist[1])));
      if(!finite(sample_moments.back())){
	sample_moments.pop_back();
	break;
      }
      indx++;
    }

    // compute sample 3 term recurrences
    MomentSequence sample_momseq(sample_moments);
    vector<double> sample_alpha;
    vector<double> sample_beta;
    sample_momseq.full_3term_recurrence(false, sample_alpha, sample_beta);


    // update average 3 term recurrences
    for(size_t i = 0; i < sample_alpha.size(); i++)
      mean_alpha[i] = mean_alpha[i]*iter/(1.0 + iter) + sample_alpha[i]/(1.0 + iter);

    for(size_t i = 0; i < sample_beta.size(); i++)
      mean_beta[i] = mean_beta[i]*iter/(1.0 + iter) + sample_beta[i]/(1.0 + iter);   

  }

  // check and remember the offdiagonals of the Jacobi matrix are the sqrt of beta 
  check_three_term_relation(mean_alpha, mean_beta);
  for(size_t i = 0; i < mean_beta.size(); i++)
    mean_beta[i] = sqrt(mean_beta[i]);

  alpha.swap(mean_alpha);
  beta.swap(mean_beta);
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
    if(d[j] == 0.0 && b_bar[j] == 0.0){
      sin_theta[j] = 0.0;
      cos_theta[j] = 1.0;
    }
    else {
      sin_theta[j] = d[j]/sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);
      cos_theta[j] = b_bar[j]/sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);
    }

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


  if(VERBOSE){
    cerr << "QR" << endl;
    cerr << "alpha = ";
    for(size_t i = 0; i < a.size(); i++)
      cerr << setprecision(16) << a[i] << ", ";
    cerr << endl;
    cerr << "beta = ";
    for(size_t i = 0; i < b.size(); i++)
      cerr << setprecision(16) << b[i] << ", ";
    cerr << endl;
  }
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
      error += fabs(qr_beta[i]);
    size_t iter = 0;
    while(iter < max_iter && error > tol){
      QRiteration(eigenvals, qr_beta, eigenvec);

      error = 0.0;
      for(size_t i = 0; i < qr_beta.size(); i++)
	error += fabs(qr_beta[i]);
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

static void
NegBin_3term_recurrence(const size_t n_points,
			const double estimated_mu,
			const double estimated_alpha,
			vector<double> &a,
			vector<double> &b){
  a.clear();
  b.clear();

  const double k = 1.0/estimated_alpha;
  const double phi = (estimated_mu + k)/estimated_mu;
  for(size_t i = 0; i < n_points; i++)
    a.push_back((2.0*i + 1.0 + k)/phi);

  for(size_t i = 1; i < n_points; i++)
    b.push_back(sqrt((i + k)*i/(phi*phi)));
}

void
MomentSequence::NegBin_quadrature_rules(const bool VERBOSE,
					const size_t n_points,
					const double tol, 
					const size_t max_iter,
					const double estimated_mu,
					const double estimated_alpha,
					vector<double> &points,
					vector<double> &weights){

  // make sure that points.size() will be less than n_points
  vector<double> a, b;
  NegBin_3term_recurrence(n_points, estimated_mu, 
			  estimated_alpha, a, b);


  if(VERBOSE){
    cerr << "QR" << endl;
    cerr << "alpha = ";
    for(size_t i = 0; i < a.size(); i++)
      cerr << setprecision(16) << a[i] << ", ";
    cerr << endl;
    cerr << "beta = ";
    for(size_t i = 0; i < b.size(); i++)
      cerr << setprecision(16) << b[i] << ", ";
    cerr << endl;
  }
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
      error += fabs(qr_beta[i]);
    size_t iter = 0;
    while(iter < max_iter && error > tol){
      QRiteration(eigenvals, qr_beta, eigenvec);

      error = 0.0;
      for(size_t i = 0; i < qr_beta.size(); i++)
	error += fabs(qr_beta[i]);
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
    }
    else{
      a.pop_back();
      b.pop_back();
    }

  }

  // now with fixed points, solve for the weights using moment equations
  // linear system in the weights
  gsl_matrix *M = gsl_matrix_alloc(points.size(), points.size());
  for(size_t i = 0; i < points.size(); i++)
    for(size_t j = 0; j < points.size(); j++)
      gsl_matrix_set(M, i, j, pow(points[j], i));
 

  int signum = 0;

  gsl_permutation *perm = gsl_permutation_alloc(points.size());

  gsl_vector *w = gsl_vector_alloc(points.size());

  gsl_vector *moms = gsl_vector_alloc(points.size());
  for(size_t i = 0; i < points.size(); i++)
    gsl_vector_set(moms, i, moments[i]);

  if(VERBOSE){
    cerr << "equation matrix:" << endl;
    for(size_t i = 0; i < points.size(); i++){
      for(size_t j = 0; j < points.size(); j++)
	cerr << gsl_matrix_get(M, i, j) << "\t";
      cerr << " = " << gsl_vector_get(moms, i) << endl;
    }
  }

  gsl_linalg_LU_decomp(M, perm, &signum);

  gsl_linalg_LU_solve(M, perm, moms, w);

  weights.resize(points.size());
  for(size_t i = 0; i < weights.size(); i++)
    weights[i] = gsl_vector_get(w, i);
}

