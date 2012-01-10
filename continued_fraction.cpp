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

#include <vector>
#include <cmath>
#include <cassert>
#include <complex>

using std::vector;
using std::complex;
using std::real;
using std::imag;

/*
// This is the quotient difference algorithm...
void
quotient_difference_alg_for_contfrac_coeffs(const vector<double> &coeffs,
const size_t depth,
vector<double> &cf_coeffs){
vector< vector<double> > q_table(depth, vector<double>(depth, 0.0));
vector< vector<double> > e_table(depth, vector<double>(depth, 0.0));
for(size_t i = 0; i < q_table[1].size(); i++)
q_table[1][i] = coeffs[i+1]/coeffs[i];
  
for(size_t j = 0; j < depth-1; j++)
e_table[1][j] = q_table[1][j+1] - q_table[1][j] + e_table[0][j+1];
  
for(size_t i = 2; i < depth; i++){
for(size_t j = 0; j < depth - i; j++)
q_table[i][j] = q_table[i-1][j+1]*e_table[i-1][j+1]/e_table[i-1][j];
    
for(size_t j = 0; j < depth - i; j++)
e_table[i][j] = q_table[i][j+1] - q_table[i][j] + e_table[i-1][j+1];
}
  
//   for (size_t j = 0; j < q_table.front().size(); ++j) {
//     for (size_t i = 0; i < q_table.size(); ++i)
//       std::cerr << std::setprecision(5) << std::setw(10) << q_table[j][i] << ' ';
//     std::cerr << std::endl;
//   }
//   std::cerr << std::endl;
//   for (size_t j = 0; j < e_table.front().size(); ++j) {
//     for (size_t i = 0; i < e_table.size(); ++i)
//       std::cerr << std::setprecision(5) << std::setw(10) << e_table[j][i] << ' ';
//     std::cerr << std::endl;
//   }
//   std::cerr << std::endl;
  
cf_coeffs.push_back(coeffs[0]);
for (size_t i = 1; i < depth; ++i) {
if (i % 2 == 0) cf_coeffs.push_back(-e_table[i/2][0]);
else cf_coeffs.push_back(-q_table[(i + 1)/2][0]);
}
}


void 
product_difference_alg_for_contfrac_coeffs(const vector<double> &coeffs, 
const size_t depth,
vector<double> &cf_coeffs) {
vector< vector<double> > p_table(2*depth-2, vector<double>(2*depth-2, 0.0));
p_table[0][0] = 1;
for(size_t i = 0; i < p_table.size(); i++)
p_table[i][1] = coeffs[i];
for(size_t j = 2; j < p_table[0].size(); j++){
for(size_t i = 0; i < p_table.size()-1; i++){
p_table[i][j] = p_table[0][j-1]*p_table[i+1][j-2] - p_table[0][j-2]*p_table[i+1][j-1];
}
}
cf_coeffs.push_back(coeffs[0]);
for(size_t i = 1; i < depth; i++)
cf_coeffs.push_back(p_table[0][i+1]/(p_table[0][i-1]*p_table[0][i]));
}
*/

// quotient-difference algorithm to compute ContFrac coeffs
static void
ContinuedFraction_qd(const vector<double> &coeffs, 
		     vector<double> &cf_coeffs){
  const size_t depth = coeffs.size();
  vector< vector<double> > q_table(depth, vector<double>(depth, 0.0));
  vector< vector<double> > e_table(depth, vector<double>(depth, 0.0));
  for(size_t i = 0; i < q_table[1].size(); i++)
    q_table[1][i] = coeffs[i + 1]/coeffs[i];
  
  for(size_t j = 0; j < depth-1; j++)
    e_table[1][j] = 
      q_table[1][j + 1] - q_table[1][j] + e_table[0][j + 1];
  
  for(size_t i = 2; i < depth; i++){
    for(size_t j = 0; j < depth; j++)
      q_table[i][j] = 
	q_table[i - 1][j + 1]*e_table[i - 1][j + 1]/e_table[i - 1][j];
    for(size_t j = 0; j < depth; j++)
      e_table[i][j] = 
	q_table[i][j + 1] - q_table[i][j] + e_table[i - 1][j + 1];
  }
  if (CF_APPROX_VERBOSE)
    cerr << "q calculated" << endl;
  cf_coeffs.push_back(coeffs[0]);
  for(size_t i = 1; i < depth; i++){
    if(i % 2 == 0)
      cf_coeffs.push_back(-e_table[i/2][0]);
    else
      cf_coeffs.push_back(-q_table[(i + 1)/2][0]);
  }
}

// compute CF coeffs when upper_offset > 0
static void
ContinuedFraction_upper_offset(const vector<double> &coeffs,
			       const size_t depth,
			       const size_t offset,
			       vector<double> &offset_cf_coeffs,
			       vector<double> &cf_coeffs){ 
//first offset coefficients set to first offset coeffs
  vector<double> holding_coeffs;
  for(size_t i = offset; i < depth; i++)
    holding_coeffs.push_back(coeffs[i]);
  ContinuedFraction_qd(holding_coeffs, depth-offset, cf_coeffs);
  for(size_t i = 0; i < offset; i++)
    offset_cf_coeffs.push_back(coeffs[i]);
}

// evaluate CF when upper_offset > 0
// using euler's recursion
static double
cf_approx_upper_offset(const vector<double> &cf_coeffs,
		       const vector<double> &offset_coeffs,
		       const double t,
		       const double tolerance){
  if(t == 0)
    return 0.0;
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*t*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*t*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      //rescale to avoid over- and underflow
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
    double offset_terms = 0.0;
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(t, i);
    return t*(offset_terms + 
	      pow(t, offset_coeffs.size())*current_num/current_denom);
  }
} 

// calculate CF coeffs when lower_offset > 0
static void
ContinuedFraction_lower_offset(const vector<double> &coeffs,
			       const size_t offset,
			       vector<double> &offset_cf_coeffs,
			       vector<double> &cf_coeffs){ 
  //need to work with reciprocal series g = 1/f, then invert
  //set offset_coeffs to 1st offset coeffs of 1/f
  vector<double> reciprocal_coeffs;
  reciprocal_coeffs.push_back(1/coeffs[0]);
  for(size_t i = 1; i < depth; i++){
    double holding_val = 0.0;
    for(size_t j = 0; j < i; j++)
      holding_val += coeffs[i - j]*reciprocal_coeffs[j];
    
    reciprocal_coeffs.push_back(-holding_val/coeffs[0]);
  }
  vector<double> holding_coeffs;
  for(size_t i = offset; i < coeffs.size(); i++)
    holding_coeffs.push_back(reciprocal_coeffs[i]);
  ContinuedFraction_qd(holding_coeffs, cf_coeffs);
  for(size_t i = 0; i < offset; i++)
    offset_cf_coeffs.push_back(reciprocal_coeffs[i]);
}

// calculate cont_frac approx when lower_offdiag > 0
static double 
cf_approx_lower_offset(const vector<double> &cf_coeffs,
		       const vector<double> &offset_coeffs,
		       const double t, const double tolerance){
  if(t == 0)
    return 0.0;
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;
    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 
    for(size_t i = 1; i < cf_coeffs.size(); i++){
      current_num = prev_num1 + cf_coeffs[i]*t*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*t*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      //rescale terms to avoid over and underflow
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
    double offset_terms = 0.0;
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(t, i);
    return t/(offset_terms + 
	      pow(t, offset_coeffs.size())*current_num/current_denom);
  }
}

// calculate cont_frac approx when there is no offset
static double
cf_approx_no_offset(const vector<double> &cf_coeffs, //uses euler's recursion
		    const double t, const double tolerance){
  if(t == 0.0){
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
      current_num = prev_num1 + cf_coeffs[i]*t*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*t*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      //rescale to avoid over and underflow
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
    return t*current_num/current_denom;
  }
}

// calculate cf_coeffs depending on offset
static void
compute_cf_coeffs(ContinuedFraction &cont_frac) {
  cf_coeffs.clear();
  offset_coeffs.clear();
  if(upper_offset == 0 && lower_offset == 0){    
    vector<double> temp_cf_coeffs;
    ContinuedFraction_qd(ps_coeffs, temp_cf_coeffs);
    cont_frac.set_cf_coeffs(temp_cf_coeffs);
  }
  else if(upper_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    ContinuedFraction_upper_offset(ps_coeffs, depth, upper_offset, 
				   temp_cf_coeffs, temp_offset_coeffs);
    cont_frac.set_offset_coeffs(temp_offset_coeffs);
    cont_frac.set_cf_coeffs(temp_cf_coeffs);
  }
  else if(lower_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    ContinuedFraction_lower_offset(ps_coeffs, depth, lower_offset,
				   temp_offset_coeffs, temp_cf_coeffs);
    cont_frac.set_offset_coeffs(temp_offset_coeffs);
    cont_frac.set_cf_coeffs(temp_cf_coeffs);
  }
}

// calculate cont_frac approx depending on offset
double
ContinuedFraction::cf_approx(const double t, const double tolerance){
  //no offset
  if (upper_offset == 0 && lower_offset == 0)
    return cf_approx_no_offset(cf_coeffs, t, tolerance);
  //upper offset
  else if (upper_offset > 0)
    return cf_approx_uooer_offset(cf_coeffs, offset_coeffs,
				  t, tolerance);
  //lower offset
  else if (lower_offset > 0)
    return(cf_approx_lower_offset(cf_coeffs, offset_coeffs,
				  t, tolerance));
}

// compute cf_approx for complex values to compute deriv when no offset
static void
cf_approx_complex_no_offset(const vector<double> &cf_coeffs,
			    const complex<double> perturbed_val,
			    const double tolerance,
			    complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    
    // Previous elements of the table to recursively fill it
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);
    
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      // rescale to avoid over and underflow
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
    approx = perturbed_val*current_num/current_denom;
  }
}

// compute complex cf_approx when upper_offset > 0
static void
cf_approx_complex_upper_offset(const vector<double> &cf_coeffs,
			       const vector<double> &offset_coeffs,
			       const complex<double> perturbed_val,
			       const double tolerance,
			       complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      //rescale to avoid over and underflow
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
    complex<double> offset_terms(0.0, 0.0);
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = perturbed_val*(offset_terms + pow(perturbed_val, offset_coeffs.size())*current_num/current_denom);
  }
} 

// compute cf approx when lower_offset > 0
static void
cf_approx_complex_lower_offset(const vector<double> &cf_coeffs,
			       const vector<double> &offset_coeffs,
			       const complex<double> perturbed_val,
			       const double tolerance,
			       complex<double> &approx) {
  const complex<double> i(0.0,1.0);
  if (norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);
    for(size_t j = 1; j < cf_coeffs.size(); j++){
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      //rescale to avoid over and underflow
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
    complex<double> offset_terms(0.0, 0.0);
    for(size_t i = 0; i < offset_coeffs.size(); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = perturbed_val/(offset_terms + pow(perturbed_val, offset_coeffs.size())*current_num/current_denom);
  }
} 

//compute cf approx for complex depending on offset
double
ContinuedFraction::cf_deriv_complex(const double val,
				    const double dx,
				    const double tolerance){
  const complex<double> i(0.0,1.0);
  complex<double> df(0.0, 0.0);
  complex<double> value(val, 0.0);
  if(upper_offset == 0 && lower_offset == 0)
    cf_approx_complex_no_offset(cf_coeffs, value + dx*i, tolerance, df);
  else if(upper_offset > 0)
    cf_approx_complex_upper_offset(cf_coeffs, offset_coeffs,
				   value + dx*i, tolerance, df);
  else if(lower_offset > 0)
    cf_approx_complex_lower_offset(cf_coeffs, offset_coeffs,
				   value + dx*i, tolerance, df);
  return(imag(df)/dx);
}

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/std::max(a, b)); //delta
}

// locate zero deriv by bisection to find local max
static double
locate_zero_cf_deriv(const double val, const double prev_val,
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

double
ContinuedFraction::locate_local_max(const double max_time,
				    const double dx,
				    const double upper_bound,
				    const double tolerance){

}
