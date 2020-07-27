/*    Copyright (C) 2013-2015
 *                  University of Southern California and
 *                  Andrew D. Smith and Timothy Daley
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
#include <algorithm>

using std::string;
using std::vector;
using std::endl;
using std::max;
using std::cerr;
using std::setprecision;
using std::swap;
using std::find_if;
using std::transform;


/////////////////////////////////////////////////////
// test Hankel moment matrix
// ensure moment sequence is positive definite
// truncate moment sequence to ensure pos def
size_t
ensure_pos_def_mom_seq(vector <double> &moments,
                       const double tolerance,
                       const bool VERBOSE) {

  const size_t min_hankel_dim = 1;
  size_t hankel_dim = 2;
  if (moments.size() < 2*hankel_dim) {
    if (VERBOSE)
      cerr << "too few moments" << endl;
    return min_hankel_dim;
  }

  bool ACCEPT_HANKEL = true;
  while (ACCEPT_HANKEL && (2*hankel_dim - 1 < moments.size())) {

    gsl_matrix *hankel_mat = gsl_matrix_alloc(hankel_dim, hankel_dim);
    for (size_t c_idx = 0; c_idx < hankel_dim; c_idx++)
      for (size_t r_idx = 0; r_idx < hankel_dim; r_idx++)
        gsl_matrix_set(hankel_mat, c_idx, r_idx, moments[c_idx + r_idx]);

    int s = 0;
    gsl_permutation *perm = gsl_permutation_alloc(hankel_dim);
    gsl_linalg_LU_decomp(hankel_mat, perm, &s);
    const double hankel_mat_det = gsl_linalg_LU_det(hankel_mat, s);

    gsl_matrix *shift_hankel_mat = gsl_matrix_alloc(hankel_dim, hankel_dim);
    for (size_t c_idx = 0; c_idx < hankel_dim; c_idx++)
      for (size_t r_idx = 0; r_idx < hankel_dim; r_idx++)
        gsl_matrix_set(shift_hankel_mat, c_idx, r_idx,
                       moments[c_idx + r_idx + 1]);

    gsl_permutation *s_perm = gsl_permutation_alloc(hankel_dim);
    gsl_linalg_LU_decomp(shift_hankel_mat, s_perm, &s);
    const double shift_hankel_mat_det = gsl_linalg_LU_det(shift_hankel_mat, s);

    if (VERBOSE) {
      cerr << "dim" << '\t'
           << "hankel_det" << '\t'
           << "shifted_hankel_det" << endl;
      cerr << hankel_dim << '\t'
           << hankel_mat_det << '\t'
           << shift_hankel_mat_det << endl;
    }

    if (hankel_mat_det > tolerance && shift_hankel_mat_det > tolerance) {
      ACCEPT_HANKEL = true;
      hankel_dim++;
    }
    else {
      ACCEPT_HANKEL = false;
      hankel_dim--;
      moments.resize(2*hankel_dim);
      return hankel_dim;
    }
  }

  return max(hankel_dim - 1, min_hankel_dim);
}


/////////////////////////////////////////////////////
// 3 term relations

// check 3 term recurrence to avoid non-positive elements
// truncate if non-positive element found
static void
check_three_term_relation(vector<double> &a, vector<double> &b) {

  // abort if first entry is zero or smaller
  if (a[0] <= 0.0) {
    a.clear();
    b.clear();
  }

  for (size_t i = 0; i < b.size(); i++)
    // ADS: some strange logic here
    if (b[i] <= 0.0 || !isfinite(b[i]) ||
        a[i + 1] <= 0.0 || !isfinite(a[i + 1])) {
      b.resize(i);
      a.resize(i + 1);
      break;
    }
}

// check the moment sequence to avoid non-positive elements and
// truncate at first non-positive element if found
static void
check_moment_sequence(vector<double> &obs_moms) {

  if (obs_moms[0] <= 0.0 || !isfinite(obs_moms[0]))
    obs_moms.clear();

  for (size_t i = 1; i < obs_moms.size(); i++) {
    if (obs_moms[i] <= 0.0 || !isfinite(obs_moms[i])) {
      obs_moms.resize(i + 1);
      break;
    }
  }
}


void
MomentSequence::unmodified_Chebyshev(const bool VERBOSE) {

  const size_t n_points = static_cast<size_t>(floor(moments.size()/2));
  vector<double> a(n_points, 0.0);
  vector<double> b(n_points - 1, 0.0);

  vector< vector<double> > sigma(2*n_points, vector<double>(2*n_points, 0.0));
  // initialization
  a[0] = moments[1]/moments[0];
  // sigma[-1][l] = 0
  for (size_t l = 0; l < 2*n_points; l++)
    sigma[0][l] = moments[l];

  for (size_t k = 1; k <= n_points; k++) {
    for (size_t l = k; l < 2*n_points - k; l++) {
      sigma[k][l] = sigma[k-1][l+1] - a[k-1]*sigma[k-1][l];
      if (k > 1)
        sigma[k][l] -= b[k-2]*sigma[k-2][l];
    }
    if (k != n_points) {
      a[k] = sigma[k][k+1]/sigma[k][k] - sigma[k-1][k]/sigma[k-1][k-1];
      b[k-1] = sigma[k][k]/sigma[k-1][k-1];
    }
  }

  alpha = a;
  beta = b;
}

// un-normalized 3 term recurrence
void
MomentSequence::full_3term_recurrence(const bool VERBOSE,
                                      vector<double> &full_alpha,
                                      vector<double> &full_beta) {

  const size_t n_points = static_cast<size_t>(floor(moments.size()/2));
  vector<double> a(n_points, 0.0);
  vector<double> b(n_points - 1, 0.0);

  vector< vector<double> > sigma(2*n_points, vector<double>(2*n_points, 0.0));
  // initialization
  a[0] = moments[1]/moments[0];
  // sigma[-1][l] = 0
  for (size_t l = 0; l < 2*n_points; l++)
    sigma[0][l] = moments[l];

  for (size_t k = 1; k <= n_points; k++) {
    for (size_t l = k; l < 2*n_points - k; l++) {
      sigma[k][l] = sigma[k-1][l+1] - a[k-1]*sigma[k-1][l];
      if (k > 1)
        sigma[k][l] -= b[k-2]*sigma[k-2][l];
    }
    if (k != n_points) {
      a[k] = sigma[k][k+1]/sigma[k][k] - sigma[k-1][k]/sigma[k-1][k-1];
      b[k-1] = sigma[k][k]/sigma[k-1][k-1];
    }
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


/////////////////////////////////////////////////////
// Quadrature Methods

// one iteration of QR:
// following eq's 3.3 of Golub & Welsh
// one iteration is Z_N-1*Z_N-2*...*Z_1*X*Z_1*...*Z_N-1
// Z_j is givens matrix to zero out the j+1,j'th element of X
static void
QRiteration(vector<double> &alpha, vector<double> &beta,
            vector<double> &weights) {

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

  for (size_t j = 0; j < alpha.size() - 1; j++) {

    // for d and b_bar, j here is j-1 in G&W
    if (d[j] == 0.0 && b_bar[j] == 0.0) {
      sin_theta[j] = 0.0;
      cos_theta[j] = 1.0;
    }
    else {
      sin_theta[j] = d[j]/sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);
      cos_theta[j] = b_bar[j]/sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);
    }

    a[j] = (a_bar[j]*cos_theta[j]*cos_theta[j] +
            2*b_tilde[j]*cos_theta[j]*sin_theta[j] +
            alpha[j+1]*sin_theta[j]*sin_theta[j]);

    a_bar[j+1] = (a_bar[j]*sin_theta[j]*sin_theta[j] -
                  2*b_tilde[j]*cos_theta[j]*sin_theta[j] +
                  alpha[j+1]*cos_theta[j]*cos_theta[j]);

    if (j != 0)
      b[j-1] = sqrt(d[j]*d[j] + b_bar[j]*b_bar[j]);

    b_bar[j+1] =
      ((a_bar[j] - alpha[j+1])*sin_theta[j]*cos_theta[j] +
       b_tilde[j]*(sin_theta[j]*sin_theta[j] - cos_theta[j]*cos_theta[j]));

    b_tilde[j+1] = -beta[j+1]*cos_theta[j];

    d[j+1] = beta[j+1]*sin_theta[j];

    z[j] = z_bar[j]*cos_theta[j] + weights[j+1]*sin_theta[j];

    z_bar[j+1] = z_bar[j]*sin_theta[j] - weights[j+1]*cos_theta[j];
  }

  // last entries set equal to final "holding" values
  a.back() = a_bar.back();
  b.back() = b_bar.back();
  z.back() = z_bar.back();

  swap(alpha, a);
  swap(beta, b);
  swap(weights, z);
}


static bool
check_positivity(const vector<double> &v) {
  return find_if(begin(v), end(v),
                 [](const double x) {return x <= 0.0 || isinf(x);}) == end(v);
}


bool
MomentSequence::Lower_quadrature_rules(const bool VERBOSE,
                                       const size_t n_points,
                                       const double tol,
                                       const size_t max_iter,
                                       vector<double> &points,
                                       vector<double> &weights) {

  // make sure that points.size() will be less than n_points
  vector<double> a(alpha);
  a.resize((n_points < alpha.size()) ? n_points : alpha.size());
  vector<double> b(beta);
  b.resize((n_points - 1 < beta.size()) ? n_points - 1 : beta.size());

  check_three_term_relation(a, b);

  // See Gautschi pgs 10-13,
  // the nu here is the square of the off-diagonal
  // of the Jacobi matrix
  for (size_t i = 0; i < b.size(); i++)
    b[i] = sqrt(b[i]);

  vector<double> eigenvec(a.size(), 0.0);
  eigenvec[0] = 1.0;
  vector<double> eigenvals(a);
  vector<double> qr_beta(b);

  // in QR, off-diagonals go to zero use off diags for convergence
  double error_sum = 0.0;
  for (size_t i = 0; i < qr_beta.size(); i++)
    error_sum += fabs(qr_beta[i]);

  size_t iter = 0;
  while (iter < max_iter && error_sum > tol) {
    QRiteration(eigenvals, qr_beta, eigenvec);

    error_sum = 0.0;
    for (size_t i = 0; i < qr_beta.size(); i++)
      error_sum += fabs(qr_beta[i]);
    iter++;
  }

  // eigenvalues are on diagonal of J
  const bool points_are_positive = check_positivity(eigenvals);
  if (points_are_positive) {
    swap(points, eigenvals);
    swap(weights, eigenvec);
  }

  // square entries in the weights vector
  transform(begin(weights), end(weights), begin(weights),
            [](const double x) {return x*x;});

  return points_are_positive;
}
