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
#include <complex>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::fabs;
using std::complex;
using std::real;
using std::imag;

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

void
cont_frac_qd(const vector<double> &coeffs, //quotient-difference
             const size_t depth,  // odd depth gives M,M pade even depth gives M-1,M pade
             vector<double> &cf_coeffs){
  vector< vector<double> > q_table(depth, vector<double>(depth, 0.0));
  vector< vector<double> > e_table(depth, vector<double>(depth, 0.0));
  for(size_t i = 0; i < q_table[1].size(); i++)
    q_table[1][i] = coeffs[i+1]/coeffs[i];
  
  for(size_t j = 0; j < depth-1; j++)
    e_table[1][j] = q_table[1][j+1] - q_table[1][j] + e_table[0][j+1];
  
  for(size_t i = 2; i < depth; i++){
    for(size_t j = 0; j < depth; j++)
      q_table[i][j] = q_table[i-1][j+1]*e_table[i-1][j+1]/e_table[i-1][j];
    
    for(size_t j = 0; j < depth; j++)
      e_table[i][j] = q_table[i][j+1] - q_table[i][j] + e_table[i-1][j+1];
  }
  if (CF_APPROX_VERBOSE)
    cerr << "q calculated" << endl;
  cf_coeffs.push_back(coeffs[0]);
  for(size_t i = 1; i < depth; i++){
    if(i % 2 == 0)
      cf_coeffs.push_back(-e_table[i/2][0]);
    else
      cf_coeffs.push_back(-q_table[(i+1)/2][0]);
  }
}

void
cont_frac_upper_offdiagonal(const vector<double> &coeffs,
                            const size_t depth,
                            const size_t offset,
                            vector<double> &offset_cf_coeffs,
                            vector<double> &cf_coeffs){ //first offset coefficients set to coeffs
  vector<double> holding_coeffs;
  for(size_t i = offset; i < depth; i++)
    holding_coeffs.push_back(coeffs[i]);
  cont_frac_qd(holding_coeffs, depth-offset, cf_coeffs);
  for(size_t i =0; i < offset; i++)
    offset_cf_coeffs.push_back(coeffs[i]);
}


static double
compute_upper_offdiag_cf_approx(const vector<double> &cf_coeffs,
                                const vector<double> &offset_coeffs,
                                const double time,
                                const double tolerance){
  if(time == 0)
    return 0.0;
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      if(i % 10 == 0){ //rescale every 10th iter
        double rescale_val = fabs(current_num) + fabs(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    double offset_terms = 0.0;
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(time, i);
    return(time*(offset_terms + pow(time, offset_coeffs.size())*current_num/current_denom));
  }
} 

void
cont_frac_lower_offdiagonal(const vector<double> &coeffs,
                            const size_t depth,
                            const size_t offset,
                            vector<double> &offset_cf_coeffs,
                            vector<double> &cf_coeffs){  //remember to invert the formulas.
  //need to work with reciprocal series g = f^{-1}
  vector<double> reciprocal_coeffs;
  reciprocal_coeffs.push_back(1/coeffs[0]);
  for(size_t i = 1; i < depth; i++){
    double holding_val = 0.0;
    for(size_t j = 0; j < i; j++)
      holding_val += coeffs[i-j]*reciprocal_coeffs[j];
    
    reciprocal_coeffs.push_back(-holding_val/coeffs[0]);
  }
  vector<double> holding_coeffs;
  for(size_t i = offset; i < depth; i++)
    holding_coeffs.push_back(reciprocal_coeffs[i]);
  cont_frac_qd(holding_coeffs, depth-offset, cf_coeffs);
  for(size_t i = 0; i < offset; i++)
    offset_cf_coeffs.push_back(reciprocal_coeffs[i]);
}

static double 
compute_lower_offdiag_cf_approx(const vector<double> &cf_coeffs,
                                const vector<double> &offset_coeffs,
                                const double time,
                                const double tolerance){
  if(time == 0)
    return 0.0;
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      if(i % 10 == 0){ //rescale every 10th iter
        double rescale_val = fabs(current_num) + fabs(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    double offset_terms = 0.0;
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(time, i);
    return(time/(offset_terms + pow(time, offset_coeffs.size())*current_num/current_denom));
  }
}

static double
compute_cf_approx_euler(const vector<double> &cf_coeffs, //uses euler's recursion
                        const double time,
                        const double tolerance){
  if(time == 0.0){
    return 0.0;
  }
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*time*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*time*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      if(i % 10 == 0){ //rescale every 10th iter
        double rescale_val = fabs(current_num) + fabs(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    return(time*current_num/current_denom);
  }
}


void
cont_frac::compute_cf_coeffs(const vector<double> ps_coeffs,
                             const size_t depth){
  cf_coeffs.clear();
  offset_coeffs.clear();
  if(upper_offset != 0 && lower_offset != 0)
    cerr << "at least one offset must be zero, reset offset.\n";
  else if(upper_offset == 0 && lower_offset == 0){    
    vector<double> temp_cf_coeffs;
    cont_frac_qd(ps_coeffs, depth, temp_cf_coeffs);
    set_cf_coeffs(temp_cf_coeffs);
  }
  else if(upper_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    cont_frac_upper_offdiagonal(ps_coeffs, depth, upper_offset, 
                                temp_cf_coeffs, temp_offset_coeffs);
    set_offset_coeffs(temp_offset_coeffs);
    set_cf_coeffs(temp_cf_coeffs);
  }
  else if(lower_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    cont_frac_lower_offdiagonal(ps_coeffs, depth, lower_offset,
                                temp_offset_coeffs, temp_cf_coeffs);
    set_offset_coeffs(temp_offset_coeffs);
    set_cf_coeffs(temp_cf_coeffs);
  }
}

double
cont_frac::cf_approx(const double time, const double tolerance){
  if(upper_offset != 0 && lower_offset != 0){
    cerr << "at least one offset must be zero, reset offset.\n";
    return(0.0);
  }
  else if(upper_offset == 0 && lower_offset == 0)
    return(compute_cf_approx_euler(cf_coeffs, time, tolerance));
  else if(upper_offset > 0)
    return(compute_upper_offdiag_cf_approx(cf_coeffs, offset_coeffs,
                                           time, tolerance));
  else if(lower_offset > 0)
    return(compute_lower_offdiag_cf_approx(cf_coeffs, offset_coeffs,
                                           time, tolerance));
  else
    return 0.0;
}

static void
cf_approx_euler_complex(const vector<double> &cf_coeffs,
                        const complex<double> perturbed_val,
                        const double tolerance,
                        complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0);
    complex<double> prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0);
    complex<double> prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      if(j % 10 == 0){ //rescale every 10th iter
        double rescale_val = norm(current_num) + norm(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    approx = perturbed_val*current_num/current_denom;
  }
}

static void
compute_upper_offdiag_cf_approx_complex(const vector<double> &cf_coeffs,
                                        const vector<double> &offset_coeffs,
                                        const complex<double> perturbed_val,
                                        const double tolerance,
                                        complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0);
    complex<double> prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0);
    complex<double> prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      if(j % 10 == 0){ //rescale every 10th iter
        double rescale_val = norm(current_num) + norm(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    complex<double> offset_terms(0.0, 0.0);
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = perturbed_val*(offset_terms + pow(perturbed_val, offset_coeffs.size())*current_num/current_denom);
  }
} 

static void
compute_lower_offdiag_cf_approx_complex(const vector<double> &cf_coeffs,
                                        const vector<double> &offset_coeffs,
                                        const complex<double> perturbed_val,
                                        const double tolerance,
                                        complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0);
    complex<double> prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0);
    complex<double> prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      if(j % 10 == 0){ //rescale every 10th iter
        double rescale_val = norm(current_num) + norm(current_denom);
        if(rescale_val > 1/tolerance)
          rescale_val = 1/rescale_val;
        else if(rescale_val < tolerance)
          rescale_val = 1/rescale_val;
        else
          rescale_val = 1.0;
        current_num = current_num*rescale_val;
        current_denom = current_denom*rescale_val;
        prev_num1 = prev_num1*rescale_val;
        prev_num2 = prev_num2*rescale_val;
        prev_denom1 = prev_denom1*rescale_val;
        prev_denom2 = prev_denom2*rescale_val;
      }
    }
    complex<double> offset_terms(0.0, 0.0);
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = perturbed_val/(offset_terms + pow(perturbed_val, offset_coeffs.size())*current_num/current_denom);
  }
} 

double
cont_frac::cf_deriv_complex(const double val,
                            const double dx,
                            const double tolerance){
  const complex<double> i(0.0,1.0);
  complex<double> df(0.0, 0.0);
  complex<double> value(val, 0.0);
  if(upper_offset != 0 && lower_offset != 0)
    cerr << "at least one offset must be zero, reset offset.\n";
  if(upper_offset == 0 && lower_offset == 0)
    cf_approx_euler_complex(cf_coeffs, value + dx*i, tolerance, df);
  else if(upper_offset > 0)
    compute_upper_offdiag_cf_approx_complex(cf_coeffs, offset_coeffs,
                                            value + dx*i, tolerance, df);
  else if(lower_offset > 0)
    compute_lower_offdiag_cf_approx_complex(cf_coeffs, offset_coeffs,
                                            value + dx*i, tolerance, df);
  return(imag(df)/dx);
}

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/std::max(a, b)); //delta
}

double
cont_frac::locate_zero_cf_deriv(const double val, const double prev_val,
                                const double dx, const double tolerance){
  double val_low = prev_val;
  double deriv_low = cf_deriv_complex(val_low, dx, tolerance);
  double val_high = val;
  double deriv_high = cf_deriv_complex(val_high, dx, tolerance);
  double val_mid = val;
  double mid_deriv = 0.0;
  double diff = std::numeric_limits<double>::max();
  double prev_deriv = std::numeric_limits<double>::max();
  while(diff > tolerance && movement(val_low, val_high) > tolerance){
    val_mid = (val_low + val_high)/2.0;
    mid_deriv = cf_deriv_complex(val_mid, dx, tolerance);
    if((mid_deriv > 0 && deriv_low < 0) ||
       (mid_deriv < 0 && deriv_low > 0))
      val_high = val_mid;
    else 
      val_low = val_mid;
    deriv_low = cf_deriv_complex(val_low, dx, tolerance);
    deriv_high = cf_deriv_complex(val_high, dx, tolerance);
    diff = fabs((prev_deriv - mid_deriv)/prev_deriv);
    prev_deriv = mid_deriv;
  }
  return(val_mid);
}