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
		      const size_t n_moments,
		      std::vector<double> &weights,
		      std::vector<double> &points){

  vector<double> moment_estimates(moments);
  moment_estimates.resize(2*n_moments);

  vector<vector<double> > matrix(n_moments, vector<double>(n_moments, 0.0));
  for (size_t i = 0; i < n_moments; ++i) {
    for (size_t j = 0; j < n_moments; ++j) {
      matrix[i][j] = moment_estimates[j + i];
    }
  }

  vector<double> c;
  solve_linear_system(matrix, moment_estimates, c);
  c.push_back(1);

   gsl_poly_complex_workspace *w =
      gsl_poly_complex_workspace_alloc(n_moments + 1);
    
   vector<double> x(2*(n_moments + 1), 0.0);
   gsl_poly_complex_solve(&c[0], c.size(), w, &x[0]);

   weights.swap(c);
   points.swap(x);
   gsl_poly_complex_workspace_free (w);
}
