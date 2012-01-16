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
using std::cerr;
using std::endl;
using std::min;

const double TOLERANCE = 1e-20;
const double DERIV_DELTA = 1e-8;
const double CF_APPROX_VERBOSE = false;

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
			       const size_t offset,
			       vector<double> &offset_cf_coeffs,
			       vector<double> &cf_coeffs){ 
//first offset coefficients set to first offset coeffs
  vector<double> holding_coeffs;
  for(size_t i = offset; i < coeffs.size(); i++)
    holding_coeffs.push_back(coeffs[i]);
  // qd to determine cf_coeffs
  ContinuedFraction_qd(holding_coeffs, cf_coeffs);
  for(size_t i = 0; i < offset; i++)
    offset_cf_coeffs.push_back(coeffs[i]);
}

// evaluate CF when upper_offset > 0
// using euler's recursion
static double
ContFrac_eval_upper_offset(const vector<double> &cf_coeffs,
		       const vector<double> &offset_coeffs,
			   const double val, const size_t depth){
  if(val == 0)
    return 0.0;
  else{
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;

    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 

    for(size_t i = 1; i < depth; i++){
      // initialize
      current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;

      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;

      //rescale to avoid over- and underflow
      double rescale_val = fabs(current_num) + fabs(current_denom);
      if(rescale_val > 1/TOLERANCE)
	rescale_val = 1/rescale_val;
      else if(rescale_val < TOLERANCE)
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
    double offset_part = 0.0;
    for(size_t i = 0; i < min(offset_coeffs.size(),
				   depth); i++)
      offset_part += offset_coeffs[i]*pow(val, i);

    return val*(offset_part + 
		pow(val, min(offset_coeffs.size(),
			     depth))*current_num/current_denom);
  }
} 

// calculate CF coeffs when lower_offset > 0
static void
ContinuedFraction_lower_offset(const vector<double> &coeffs,
			       const size_t offset,
			       vector<double> &offset_cf_coeffs,
			       vector<double> &cf_coeffs){ 
  //need to work with reciprocal series g = 1/f, then invert
  vector<double> reciprocal_coeffs;
  reciprocal_coeffs.push_back(1/coeffs[0]);
  for(size_t i = 1; i < coeffs.size(); i++){
    double holding_val = 0.0;
    for(size_t j = 0; j < i; j++)
      holding_val += coeffs[i - j]*reciprocal_coeffs[j];
    
    reciprocal_coeffs.push_back(-holding_val/coeffs[0]);
  }

  //set offset_coeffs to 1st offset coeffs of 1/f 
  for(size_t i = 0; i < offset; i++)
    offset_cf_coeffs.push_back(reciprocal_coeffs[i]);
  // qd to compute cf_coeffs using remaining coeffs
  vector<double> holding_coeffs;
  for(size_t i = offset; i < coeffs.size(); i++)
    holding_coeffs.push_back(reciprocal_coeffs[i]);
  ContinuedFraction_qd(holding_coeffs, cf_coeffs);
}

// calculate cont_frac approx when lower_offdiag > 0
static double 
ContFrac_eval_lower_offset(const vector<double> &cf_coeffs,
		       const vector<double> &offset_coeffs,
			   const double val, const size_t depth) {
  if(val == 0)
    return 0.0;
  else{
    //initialize
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;

    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 

    for(size_t i = 1; i < depth; i++){
      // recursion
      current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;

      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      //rescale terms to avoid over and underflow
      double rescale_val = fabs(current_num) + fabs(current_denom);
      if(rescale_val > 1/TOLERANCE)
	rescale_val = 1/rescale_val;
      else if(rescale_val < TOLERANCE)
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
    for(size_t i = 0; i < min(offset_coeffs.size(),
				   depth); i++)
      offset_terms += offset_coeffs[i]*pow(val, i);
    // recall that if lower_offset > 0, we are working with 1/f, invert approx
    return val/(offset_terms + 
		pow(val, min(offset_coeffs.size(),
			     depth))*current_num/current_denom);
  }
}

// calculate cont_frac approx when there is no offset
// uses euler's recursion
static double
ContFrac_eval_no_offset(const vector<double> &cf_coeffs, 
			const double val, const size_t depth){
  if(val == 0.0){
    return 0.0;
  }
  else{
    // initialize
    double current_num = 0.0;
    double prev_num1 = cf_coeffs[0];
    double prev_num2 = 0.0;

    double current_denom = 0.0;
    double prev_denom1 = 1.0;
    double prev_denom2 = 1.0; 

    for(size_t i = 1; i < depth; i++){
      // recursion
      current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
      current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;

      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2= prev_denom1;
      prev_denom1 = current_denom;
      //rescale to avoid over and underflow
      double rescale_val = fabs(current_num) + fabs(current_denom);
      if(rescale_val > 1/TOLERANCE)
	rescale_val = 1/rescale_val;
      else if(rescale_val < TOLERANCE)
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
    return val*current_num/current_denom;
  }
}

// calculate cf_coeffs depending on offset
void
ContFracApprox::compute_cf_coeffs() {
  cont_frac_estimate.cf_coeffs.clear();
  cont_frac_estimate.offset_coeffs.clear();
  if(cont_frac_estimate.upper_offset == 0 && 
     cont_frac_estimate.lower_offset == 0){    
    vector<double> temp_cf_coeffs;
    ContinuedFraction_qd(cont_frac_estimate.ps_coeffs, temp_cf_coeffs);
    cont_frac_estimate.cf_coeffs = temp_cf_coeffs;
  }
  else if(cont_frac_estimate.upper_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    ContinuedFraction_upper_offset(cont_frac_estimate.ps_coeffs, 
				   cont_frac_estimate.upper_offset, 
				   temp_cf_coeffs, temp_offset_coeffs);
    cont_frac_estimate.offset_coeffs = temp_offset_coeffs;
    cont_frac_estimate.cf_coeffs = temp_cf_coeffs;
  }
  else if(cont_frac_estimate.lower_offset > 0){
    vector<double> temp_cf_coeffs;
    vector<double> temp_offset_coeffs;
    ContinuedFraction_lower_offset(cont_frac_estimate.ps_coeffs, 
				   cont_frac_estimate.lower_offset,
				   temp_offset_coeffs, temp_cf_coeffs);
    cont_frac_estimate.offset_coeffs = temp_offset_coeffs;
    cont_frac_estimate.cf_coeffs = temp_cf_coeffs;
  }
}

// calculate cont_frac approx depending on offset
double
cont_frac::evaluate(const double val, const size_t depth){
  //upper offset
  if (upper_offset > 0)
    return ContFrac_eval_upper_offset(cf_coeffs, offset_coeffs, 
				      val, depth);
  //lower offset
  else if (lower_offset > 0)
    return ContFrac_eval_lower_offset(cf_coeffs, offset_coeffs, 
				      val, depth);
  // no offset
  else
    return ContFrac_eval_no_offset(cf_coeffs, val, depth);

}

// compute ContFrac_eval for complex values to compute deriv when no offset
static void
ContFrac_eval_complex_no_offset(const vector<double> &cf_coeffs,
				const complex<double> perturbed_val,
				const size_t depth,
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
    
    for(size_t j = 1; j < depth; j++){
      //euler's recursion
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;
      // rescale to avoid over and underflow
      double rescale_val = norm(current_num) + norm(current_denom);
      if(rescale_val > 1/TOLERANCE)
	rescale_val = 1/rescale_val;
      else if(rescale_val < TOLERANCE)
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

// compute complex ContFrac_eval when upper_offset > 0
static void
ContFrac_eval_complex_upper_offset(const vector<double> &cf_coeffs,
				   const vector<double> &offset_coeffs,
				   const complex<double> perturbed_val,
				   const size_t depth,
				   complex<double> &approx){
  const complex<double> i(0.0,1.0);
  if(norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    //initialize
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);

    for(size_t j = 1; j < depth; j++){
      //eulers recursion
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;

      //rescale to avoid over and underflow
      double rescale_val = norm(current_num) + norm(current_denom);
      if(rescale_val > 1/TOLERANCE)
	rescale_val = 1/rescale_val;
      else if(rescale_val < TOLERANCE)
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
    for(size_t i = 0; i < min(offset_coeffs.size(),
				   depth); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = 
      perturbed_val*(offset_terms 
		     + pow(perturbed_val, 
			   min(offset_coeffs.size(),
			       depth))*current_num/current_denom);
  }
} 

// compute cf approx when lower_offset > 0
static void
ContFrac_eval_complex_lower_offset(const vector<double> &cf_coeffs,
				   const vector<double> &offset_coeffs,
				   const complex<double> perturbed_val,
				   const size_t depth,
				   complex<double> &approx) {
  const complex<double> i(0.0,1.0);
  if (norm(perturbed_val) == 0.0)
    approx = 0.0*i;
  else{
    // initialize
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);

    for(size_t j = 1; j < cf_coeffs.size(); j++){
      // euler's recursion
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;

      //rescale to avoid over and underflow
      double rescale_val = norm(current_num) + norm(current_denom);
      if(rescale_val > 1/TOLERANCE)
	rescale_val = 1/rescale_val;
      else if(rescale_val < TOLERANCE)
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
    for(size_t i = 0; i < min(offset_coeffs.size(),
			      depth); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    approx = 
      perturbed_val/(offset_terms 
		     + pow(perturbed_val, 
			   min(offset_coeffs.size(),
			       depth))*current_num/current_denom);
  }
} 

//compute cf approx for complex depending on offset
// df/dx = lim_{delta -> 0} Imag(f(val+i*delta))/delta
double
cont_frac::complex_deriv(const double val, const size_t depth){
  vector<double> ContFracCoeffs;
  get_cf_coeffs(ContFracCoeffs);
  vector<double> ContFracOffCoeffs;
  get_offset_coeffs(ContFracOffCoeffs);

  const complex<double> i(0.0,1.0);
  complex<double> df(0.0, 0.0);
  complex<double> value(val, 0.0);

  if(upper_offset == 0 && lower_offset == 0)
    ContFrac_eval_complex_no_offset(ContFracCoeffs, value + DERIV_DELTA*i, 
				    depth, df);

  else if(upper_offset > 0)
    ContFrac_eval_complex_upper_offset(ContFracCoeffs, ContFracOffCoeffs,
				       value + DERIV_DELTA*i, 
				       depth, df);

  else if(lower_offset > 0)
    ContFrac_eval_complex_lower_offset(ContFracCoeffs, ContFracOffCoeffs,
				       value + DERIV_DELTA*i, 
				       depth, df);

  return(imag(df)/DERIV_DELTA);
}

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/std::max(a, b)); //delta
}

// locate zero deriv by bisection to find local max
// within (prev_val, val)
double
ContFracApprox::locate_zero_cf_deriv(const double val, 
				     const double prev_val){
  double val_low = prev_val;
  double deriv_low = complex_deriv(val_low);
  double val_high = val;
  double deriv_high = complex_deriv(val_high);
  double val_mid = (val-prev_val)/2;
  double deriv_mid = std::numeric_limits<double>::max();

  double diff = std::numeric_limits<double>::max();
  double prev_deriv = std::numeric_limits<double>::max();

  while(diff > TOLERANCE && movement(val_low, val_high) > TOLERANCE){
    val_mid = (val_low + val_high)/2.0;
    deriv_mid = complex_deriv(val_mid);

    if((deriv_mid > 0 && deriv_low < 0) ||
       (deriv_mid < 0 && deriv_low > 0))
      val_high = val_mid;
    else 
      val_low = val_mid;

    deriv_low = complex_deriv(val_low);
    deriv_high = complex_deriv(val_high);
    diff = fabs((prev_deriv - deriv_mid)/prev_deriv);
    prev_deriv = deriv_mid;
  }

  return(val_mid);
}

// 
static bool
test_stability_local_max(const double approx, const double deriv_val,
			 const double approx_upper_bound, const double deriv_upper_bound) {
  return (deriv_val < deriv_upper_bound) &&
    (approx < approx_upper_bound);
}

// search (min_val, max_val) for local max
// return location of local max
double
ContFracApprox::locate_local_max(const double min_val, const double max_val,
				 const double step_size, const double upper_bound,
				 const double deriv_upper_bound) {
  double val = min_val;
  double prev_approx = evaluate(val);
  double prev_deriv = complex_deriv(val);
  double current_approx, current_deriv;
  
  double current_max = prev_approx;
  double current_max_loc = val;

  while(val <= max_val){
    val += step_size;
    current_approx = evaluate(val);
    current_deriv = complex_deriv(val);

    // test stability to locate possible defects
    // do not use approx if estimate is not stable
    if(test_stability_local_max(current_approx, current_deriv,
				upper_bound, deriv_upper_bound)){
      // update max if it is greater
      if((current_deriv < 0.0) && (prev_deriv > 0.0)){ 
	double possible_max_loc = 
	  locate_zero_cf_deriv(val, val-step_size);
	if(evaluate(possible_max_loc) > current_max){
	  current_max = evaluate(possible_max_loc);
	  current_max_loc = possible_max_loc;
	}
      }
    }
    //exit while loop if Approx is unstable
    else{
      val = max_val; 
      current_max = evaluate(min_val);
      current_max_loc = min_val;
    }
  }

  return current_max_loc;
}

void
ContFracApprox::set_depth(const size_t max_terms){
  depth = max_terms;
  assert(depth >= MINIMUM_ALLOWED_DEGREE);
}
