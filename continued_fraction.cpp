/*    Copyright (C) 2011 University of Southern California and
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

#include "continued_fraction.hpp"

#include <smithlab_utils.hpp>

#include <gsl/gsl_vector_double.h>  
#include <gsl/gsl_matrix_double.h>  
#include <gsl/gsl_linalg.h>  

#include <iostream>
#include <cassert>
#include <numeric>
#include <vector>
#include <cstdlib>
#include <cmath>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::fabs;

bool CF_APPROX_VERBOSE = false;

void 
cont_frac_pd(const vector<double> &coeffs, const size_t depth, 
	     vector<double> &cf_coeffs) {

  vector< vector<double> > p_vec(2*depth-2, vector<double>(2*depth-2, 0.0));
  p_vec[0][0] = 1;

  for (size_t i = 0; i < p_vec.size(); i++)
    p_vec[i][1] = coeffs[i];
  
  for (size_t j = 2; j < p_vec[0].size(); j++)
    for (size_t i = 0; i < p_vec.size()-1; i++)
      p_vec[i][j] = p_vec[0][j-1]*p_vec[i+1][j-2] - p_vec[0][j-2]*p_vec[i+1][j-1];

  if (CF_APPROX_VERBOSE)
    cerr << "p calculated" << endl;
  cf_coeffs.push_back(coeffs[0]);
  for (size_t i = 1; i < depth; i++)
    cf_coeffs.push_back(p_vec[0][i + 1]/(p_vec[0][i - 1]*p_vec[0][i]));
}

double 
compute_cf_approx(const vector<double> &cf_coeffs, const double time) {
  double val = 0.0;
  for (size_t i = cf_coeffs.size(); i > 0; --i)
    val = cf_coeffs[i-1]*time/(1.0 + val);
  return val;
}

double
compute_cf_approx_euler(const vector<double> &cf_coeffs, const double time) {
  double curr_num = 0.0, curr_denom = 0.0;
  double prev_num1 = 0.0, prev_num2 = 1.0;
  double prev_denom1 = 0.0, prev_denom2 = 1.0; 
  
  for (size_t i = 0; i < cf_coeffs.size(); ++i) {
    curr_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
    curr_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
    prev_num2 = prev_num1;
    prev_num1 = curr_num;
    prev_denom2= prev_denom1;
    prev_denom1 = curr_denom;
    if (CF_APPROX_VERBOSE)
      cerr << curr_num << '\t' << curr_denom << endl;
  }
  return curr_num/curr_denom;
}
