/* Copyright (C) 2013-2026 Andrew D. Smith and Timothy Daley
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "moment_sequence.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <iterator>
#include <numeric>
#include <ranges>
#include <utility>  // IWYU pragma: keep
#include <vector>

[[nodiscard]] constexpr auto
is_positive(const double x) {
  return x > 0.0 && std::isfinite(x);
}

// check the moment sequence to avoid non-positive elements and truncate at
// first non-positive element if found
static auto
truncate_before_non_positive(std::vector<double> &moments) {
  const auto lam = [&](const auto x) { return x <= 0.0 || !std::isfinite(x); };
  const auto first_non_positive = std::ranges::find_if(moments, lam);
  moments.resize(std::distance(std::begin(moments), first_non_positive));
}

MomentSequence::MomentSequence(const std::vector<double> &m) : moments{m} {
  truncate_before_non_positive(moments);
  unmodified_chebyshev();  // calculate 3-term recurrence
}

void
LU_decomp(std::vector<std::vector<double>> &A, std::vector<int> &P) {
  const std::size_t N = std::size(A);

  P.resize(N + 1);
  std::ranges::generate_n(std::begin(P), std::ssize(P),
                          [n{0}]  // cppcheck-suppress[syntaxError]
                          mutable { return n++; });

  for (auto i = 0LU; i < N; ++i) {
    double maxA{};
    std::size_t max_i = i;

    for (auto k = i; k < N; ++k) {
      auto absA = std::fabs(A[k][i]);
      if (absA > maxA) {
        maxA = absA;
        max_i = k;
      }
    }

    if (max_i != i) {
      std::swap(P[i], P[max_i]);  // pivoting P
      std::swap(A[i], A[max_i]);  // pivoting rows of A
      ++P[N];  // counting pivots starting from N (for determinant)
    }

    for (auto j = i + 1; j < N; ++j) {
      A[j][i] /= A[i][i];
      for (auto k = i + 1; k < N; ++k)
        A[j][k] -= A[j][i] * A[i][k];
    }
  }
}

[[nodiscard]] static auto
LU_determinant(const std::vector<std::vector<double>> &A,
               const std::vector<int> &P) -> double {
  assert(std::size(A) < std::size(P));
  const auto N = std::size(A);
  double det = 1.0;
  for (auto i = 0LU; i < N; ++i)
    det *= A[i][i];
  return (((P[N] - N) % 2 == 0) ? 1 : -1) * det;
}

/////////////////////////////////////////////////////
// test Hankel moment matrix

// ensure moment sequence is positive definite by truncating moment sequence as
// needed
auto
ensure_positive_definite_moment_sequence(
  std::vector<double> &moments, const double tolerance) -> std::size_t {
  const std::size_t min_hankel_dim = 1;
  std::size_t hankel_dim = 2;
  if (std::size(moments) < 2 * hankel_dim)
    return min_hankel_dim;

  while (2 * hankel_dim - 1 < std::size(moments)) {
    std::vector<std::vector<double>> hankel_mat(
      hankel_dim, std::vector<double>(hankel_dim, 0.0));
    for (std::size_t c_idx = 0; c_idx < hankel_dim; ++c_idx)
      for (std::size_t r_idx = 0; r_idx < hankel_dim; ++r_idx)
        hankel_mat[c_idx][r_idx] = moments[c_idx + r_idx];

    std::vector<int> perm;
    LU_decomp(hankel_mat, perm);
    const double hankel_mat_det = LU_determinant(hankel_mat, perm);

    std::vector<std::vector<double>> shift_hankel_matrix(
      hankel_dim, std::vector<double>(hankel_dim, 0.0));
    for (std::size_t c_idx = 0; c_idx < hankel_dim; ++c_idx)
      for (std::size_t r_idx = 0; r_idx < hankel_dim; ++r_idx)
        shift_hankel_matrix[c_idx][r_idx] = moments[c_idx + r_idx + 1];

    std::vector<int> s_perm;
    LU_decomp(shift_hankel_matrix, s_perm);
    const double shift_hankel_mat_det =
      LU_determinant(shift_hankel_matrix, s_perm);

    if (hankel_mat_det > tolerance && shift_hankel_mat_det > tolerance) {
      ++hankel_dim;
    }
    else {
      --hankel_dim;
      moments.resize(2 * hankel_dim);
      return hankel_dim;
    }
  }

  return std::max(hankel_dim - 1, min_hankel_dim);
}

/// 3 term relations

// check 3 term recurrence to avoid non-positive elements truncate if
// non-positive element found
static void
check_three_term_relation(std::vector<double> &a, std::vector<double> &b) {
  // abort if first entry is zero or smaller
  if (a[0] <= 0.0) {
    a.clear();
    b.clear();
  }
  for (std::size_t i = 0; i < std::size(b); ++i)
    // ADS: some strange logic here
    if (b[i] <= 0.0 || !std::isfinite(b[i]) || a[i + 1] <= 0.0 ||
        !std::isfinite(a[i + 1])) {
      b.resize(i);
      a.resize(i + 1);
      break;
    }
}

void
MomentSequence::unmodified_chebyshev() {
  const auto n_points = static_cast<std::size_t>(
    std::floor(static_cast<double>(std::size(moments)) / 2.0));
  std::vector<double> a(n_points, 0.0);
  std::vector<double> b(n_points - 1, 0.0);

  const auto dim = 2 * n_points;
  auto sigma = std::vector(dim, std::vector<double>(dim, 0.0));
  // initialization
  a[0] = moments[1] / moments[0];
  // sigma[-1][l] = 0
  for (std::size_t l = 0; l < dim; ++l)
    sigma[0][l] = moments[l];

  for (std::size_t k = 1; k <= n_points; ++k) {
    for (std::size_t l = k; l < dim - k; ++l) {
      sigma[k][l] = sigma[k - 1][l + 1] - a[k - 1] * sigma[k - 1][l];
      if (k > 1)
        sigma[k][l] -= b[k - 2] * sigma[k - 2][l];
    }
    if (k != n_points) {
      a[k] =
        sigma[k][k + 1] / sigma[k][k] - sigma[k - 1][k] / sigma[k - 1][k - 1];
      b[k - 1] = sigma[k][k] / sigma[k - 1][k - 1];
    }
  }

  alpha = a;
  beta = b;
}

// un-normalized 3 term recurrence
void
MomentSequence::full_three_term_recurrence(std::vector<double> &full_alpha,
                                           std::vector<double> &full_beta) {
  const auto n_points = static_cast<std::size_t>(
    std::floor(static_cast<double>(std::size(moments)) / 2.0));

  std::vector<double> a(n_points, 0.0);
  std::vector<double> b(n_points - 1, 0.0);

  const auto dim = 2 * n_points;
  auto sigma = std::vector(dim, std::vector<double>(dim, 0.0));
  // initialization
  a[0] = moments[1] / moments[0];
  // sigma[-1][l] = 0
  for (std::size_t l = 0; l < dim; ++l)
    sigma[0][l] = moments[l];

  for (std::size_t k = 1; k <= n_points; ++k) {
    for (std::size_t l = k; l < dim - k; ++l) {
      sigma[k][l] = sigma[k - 1][l + 1] - a[k - 1] * sigma[k - 1][l];
      if (k > 1)
        sigma[k][l] -= b[k - 2] * sigma[k - 2][l];
    }
    if (k != n_points) {
      a[k] =
        sigma[k][k + 1] / sigma[k][k] - sigma[k - 1][k] / sigma[k - 1][k - 1];
      b[k - 1] = sigma[k][k] / sigma[k - 1][k - 1];
    }
  }

  full_alpha.swap(a);
  full_beta.swap(b);
}

/////////////////////////////////////////////////////
// Quadrature Methods

// one iteration of QR:
// following eq's 3.3 of Golub & Welsh
// one iteration is Z_N-1*Z_N-2*...*Z_1*X*Z_1*...*Z_N-1
// Z_j is givens matrix to zero out the j+1,j'th element of X
static void
QR_iteration(std::vector<double> &alpha, std::vector<double> &beta,
             std::vector<double> &weights) {
  // initialize variables
  std::vector<double> sin_theta(std::size(alpha), 0.0);
  std::vector<double> cos_theta(std::size(alpha), 0.0);

  std::vector<double> a(std::size(alpha), 0.0);
  std::vector<double> a_bar(std::size(alpha), 0.0);
  a_bar[0] = alpha[0];

  std::vector<double> b(beta);
  std::vector<double> b_bar(std::size(alpha), 0.0);
  b_bar[0] = alpha[0];
  std::vector<double> b_tilde(std::size(alpha), 0.0);
  b_tilde[0] = beta[0];

  std::vector<double> d(std::size(alpha), 0.0);
  d[0] = beta[0];

  std::vector<double> z(weights);
  std::vector<double> z_bar(std::size(weights), 0.0);
  z_bar[0] = z[0];

  for (std::size_t j = 0; j + 1 < std::size(alpha); ++j) {
    // for d and b_bar, j here is j-1 in G&W
    if (d[j] == 0.0 && b_bar[j] == 0.0) {
      sin_theta[j] = 0.0;
      cos_theta[j] = 1.0;
    }
    else {
      sin_theta[j] = d[j] / std::sqrt(d[j] * d[j] + b_bar[j] * b_bar[j]);
      cos_theta[j] = b_bar[j] / std::sqrt(d[j] * d[j] + b_bar[j] * b_bar[j]);
    }

    a[j] = (a_bar[j] * cos_theta[j] * cos_theta[j] +
            2 * b_tilde[j] * cos_theta[j] * sin_theta[j] +
            alpha[j + 1] * sin_theta[j] * sin_theta[j]);

    a_bar[j + 1] = (a_bar[j] * sin_theta[j] * sin_theta[j] -
                    2 * b_tilde[j] * cos_theta[j] * sin_theta[j] +
                    alpha[j + 1] * cos_theta[j] * cos_theta[j]);

    if (j != 0)
      b[j - 1] = std::sqrt(d[j] * d[j] + b_bar[j] * b_bar[j]);

    b_bar[j + 1] = ((a_bar[j] - alpha[j + 1]) * sin_theta[j] * cos_theta[j] +
                    b_tilde[j] * (sin_theta[j] * sin_theta[j] -
                                  cos_theta[j] * cos_theta[j]));

    b_tilde[j + 1] = -beta[j + 1] * cos_theta[j];

    d[j + 1] = beta[j + 1] * sin_theta[j];

    z[j] = z_bar[j] * cos_theta[j] + weights[j + 1] * sin_theta[j];

    z_bar[j + 1] = z_bar[j] * sin_theta[j] - weights[j + 1] * cos_theta[j];
  }

  // last entries set equal to final "holding" values
  a.back() = a_bar.back();
  b.back() = b_bar.back();
  z.back() = z_bar.back();

  std::swap(alpha, a);
  std::swap(beta, b);
  std::swap(weights, z);
}

[[nodiscard]] static auto
check_positivity(const std::vector<double> &v) -> bool {
  const auto not_positive = [](const auto x) {
    return x <= 0.0 || std::isinf(x);
  };
  return std::ranges::find_if(v, not_positive) == std::cend(v);
}

auto
MomentSequence::lower_quadrature_rules(const std::size_t n_points,
                                       const double tol,
                                       const std::size_t max_iter)
  -> std::tuple<std::vector<double>, std::vector<double>> {
  const auto abs_sum = [](const auto a, const auto b) {
    return a + std::fabs(b);
  };
  const auto abs_sum_vec = [&](const auto &v) {
    return std::reduce(std::cbegin(v), std::cend(v), 0.0, abs_sum);
  };
  // make sure that the number of points will be at most n_points
  auto a = alpha;
  if (n_points < std::size(a))
    a.resize(n_points);
  auto b = beta;
  if (n_points < std::size(b) + 1)
    b.resize(n_points - 1);

  check_three_term_relation(a, b);

  // See Gautschi pgs 10-13,
  // the nu here is the square of the off-diagonal of the Jacobi matrix
  std::ranges::for_each(b, [](auto &b_val) { b_val = std::sqrt(b_val); });

  std::vector<double> eigenvec(std::size(a), 0.0);
  eigenvec[0] = 1.0;
  auto eigenvals = a;
  auto qr_beta = b;

  // in QR, off-diagonals go to zero use off diags for convergence
  auto error_sum = abs_sum_vec(qr_beta);
  for (std::size_t iter{}; iter < max_iter && error_sum > tol; ++iter) {
    QR_iteration(eigenvals, qr_beta, eigenvec);
    error_sum = abs_sum_vec(qr_beta);
  }

  // eigenvalues are on diagonal of J
  const bool points_are_positive = check_positivity(eigenvals);
  if (!points_are_positive) {
    eigenvals.clear();
    eigenvec.clear();
  }

  // square entries in the weights vector
  std::ranges::for_each(eigenvec, [](auto &x) { x *= x; });

  return std::tuple{std::move(eigenvals), std::move(eigenvec)};
}
