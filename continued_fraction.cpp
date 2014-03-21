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

#include "continued_fraction.hpp"
#include <smithlab_utils.hpp>
#include <RNG.hpp>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <vector>
#include <cmath>
#include <cassert>
#include <complex>
#include <limits>
#include <sys/types.h>
#include <unistd.h>

using std::vector;
using std::complex;
using std::real;
using std::imag;
using std::cerr;
using std::endl;
using std::min;

const double TOLERANCE = 1e-20;
const double DERIV_DELTA = 1e-8;




static double
get_rescale_value(const double numerator, const double denominator) {
  const double rescale_val = fabs(numerator) + fabs(denominator);
  if (rescale_val > 1.0/TOLERANCE)
    return 1.0/rescale_val;
  else if (rescale_val < TOLERANCE)
    return 1.0/rescale_val;
  return 1.0;
}

static double
get_rescale_value(const complex<double> numerator, const complex<double> denominator) {
  const double rescale_val = norm(numerator) + norm(denominator);
  if (rescale_val > 1.0/TOLERANCE)
    return 1.0/rescale_val;
  if (rescale_val < TOLERANCE)
    return 1.0/rescale_val;
  return 1.0;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  QUOTIENT DIFFERENCE ALGORITHMS
////

/* quotient-difference algorithm to compute continued fraction
   coefficients
*/ 
static void
quotdiff_algorithm(const vector<double> &ps_coeffs, vector<double> &cf_coeffs) { //vector for power series coefficients & vector for continued fraction coefficients
  
  const size_t depth = ps_coeffs.size(); //degree of power series
  vector< vector<double> > q_table(depth, vector<double>(depth+1, 0.0));
  vector< vector<double> > e_table(depth, vector<double>(depth+1, 0.0));

  for (size_t j = 0; j < depth-1; j++) 
    q_table[1][j] = ps_coeffs[j + 1]/ps_coeffs[j];
  
  for (size_t j = 0; j < depth-1; j++) 
    e_table[1][j] = q_table[1][j + 1] - q_table[1][j] + e_table[0][j + 1];
    
  //using intial values of E(i)(j)'s and Q(i)(j)'s, fill rest of the q table and e table
  for (size_t i = 2; i < depth; i++) {
    for (size_t j = 0; j < depth; j++)
      q_table[i][j] = q_table[i - 1][j + 1]*e_table[i - 1][j + 1]/e_table[i - 1][j];
    
    for (size_t j = 0; j < depth; j++)
      e_table[i][j] = q_table[i][j + 1] - q_table[i][j] + e_table[i - 1][j + 1];
  }

  cf_coeffs.push_back(ps_coeffs[0]); //add first power series coefficient to end of vector for continued fraction coefficients

  //setting coefficients for continued fraction 
  for (size_t i = 1; i < depth; ++i) {
    if (i % 2 == 0) 
      cf_coeffs.push_back(-e_table[i/2][0]);
    else
      cf_coeffs.push_back(-q_table[(i + 1)/2][0]);
  }
}  


// compute CF coeffs when upper_offset > 0
//above the diagonal referring to degree of polynomial in numerator of Pade approximant is greater than degree of polynomial in the denominator 
static void
quotdiff_above_diagonal(const vector<double> &coeffs, const size_t offset,
                        vector<double> &cf_coeffs, vector<double> &offset_coeffs) {  
  //first offset coefficients set to first offset coeffs
  vector<double> holding_coeffs; 
  for (size_t i = offset; i < coeffs.size(); i++)
    holding_coeffs.push_back(coeffs[i]);
  
  // qd to determine cf_coeffs
  quotdiff_algorithm(holding_coeffs, cf_coeffs);
  for (size_t i = 0; i < offset; i++)
    offset_coeffs.push_back(coeffs[i]);
}


// calculate CF coeffs when lower_offset > 0
static void
quotdiff_below_diagonal(const vector<double> &coeffs, const size_t offset, 
                        vector<double> &cf_coeffs, vector<double> &offset_coeffs) {
  
  //need to work with reciprocal series g = 1/f, then invert
  vector<double> reciprocal_coeffs;
  reciprocal_coeffs.push_back(1.0/coeffs[0]);
  for (size_t i = 1; i < coeffs.size(); i++) {
    double holding_val = 0.0;
    for (size_t j = 0; j < i; ++j)
      holding_val += coeffs[i - j]*reciprocal_coeffs[j];
    reciprocal_coeffs.push_back(-holding_val/coeffs[0]);
  }
  
  //set offset_coeffs to 1st offset coeffs of 1/f 
  for (size_t i = 0; i < offset; i++)
    offset_coeffs.push_back(reciprocal_coeffs[i]);
  
  // qd to compute cf_coeffs using remaining coeffs
  vector<double> holding_coeffs;
  for (size_t i = offset; i < coeffs.size(); i++)
    holding_coeffs.push_back(reciprocal_coeffs[i]);
  
  quotdiff_algorithm(holding_coeffs, cf_coeffs);
}

// output new ContinuedFraction with decreased degree
// and coeffs equal to the old, but decreased in degree
ContinuedFraction
ContinuedFraction::decrease_degree(const ContinuedFraction &CF,
                                   const size_t decrement) {
  // create return ContinuedFraction
  ContinuedFraction decreasedCF;
  // properties of orig CF to decrement
  vector<double> decreased_ps_coeffs(CF.ps_coeffs);
  vector<double> decreased_cf_coeffs(CF.cf_coeffs);
  // decrease order
  for(size_t i = 0; i < decrement; i++) {
    decreased_ps_coeffs.pop_back();
    decreased_cf_coeffs.pop_back();
  }

  // just a copy
  vector<double> decreased_offset_coeffs(CF.offset_coeffs);

  // set return ContinuedFraction
  decreasedCF.ps_coeffs = decreased_ps_coeffs;
  decreasedCF.cf_coeffs = decreased_cf_coeffs;
  decreasedCF.offset_coeffs = decreased_offset_coeffs;
  decreasedCF.diagonal_idx = CF.diagonal_idx;
  decreasedCF.degree = CF.degree - decrement;

  return decreasedCF;
}

ContinuedFraction
ContinuedFraction::truncate_degree(const ContinuedFraction &CF,
				   const size_t n_terms){
  ContinuedFraction truncated_CF;
  if(CF.degree < n_terms){
    cerr << "current CF degree   = " << CF.degree << endl;
    cerr << "truncated CF degree = " << n_terms << endl; 
    throw SMITHLABException("degree of truncate CF must be at least as large as current");
  }

  vector<double> truncated_ps_coeffs(CF.ps_coeffs);
  vector<double> truncated_cf_coeffs(CF.cf_coeffs);
  vector<double> truncated_offset_coeffs(CF.offset_coeffs);

  truncated_ps_coeffs.resize(n_terms);
  truncated_cf_coeffs.resize(n_terms - truncated_offset_coeffs.size());

  truncated_CF.ps_coeffs = truncated_ps_coeffs;
  truncated_CF.cf_coeffs = truncated_cf_coeffs;
  truncated_CF.offset_coeffs = truncated_offset_coeffs;
  truncated_CF.diagonal_idx = CF.diagonal_idx;
  truncated_CF.degree = n_terms;

  return truncated_CF;
}

ContinuedFraction::ContinuedFraction(const vector<double> &ps_cf, 
                                     const int di, const size_t dg) :
  ps_coeffs(ps_cf), diagonal_idx(di), degree(dg) {

  if (diagonal_idx == 0)
    quotdiff_algorithm(ps_coeffs, cf_coeffs);
  else if (diagonal_idx > 0)
    quotdiff_above_diagonal(ps_coeffs, diagonal_idx, cf_coeffs, offset_coeffs);
  else // if(cont_frac_estimate.lower_offset > 0) {
    quotdiff_below_diagonal(ps_coeffs, -diagonal_idx, cf_coeffs, offset_coeffs);
  // notice the "-" above so that -diagonal_idx > 0
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  FUNCTIONS TO EVALUATE CONTINUED FRACTIONS AT A POINT
////

/* evaluate CF when upper_offset > 0 using Euler's recursion
 */
static double
evaluate_above_diagonal(const vector<double> &cf_coeffs,
                        const vector<double> &offset_coeffs,
                        const double val, const size_t depth) {
  
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;
  
  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0; 
  
  for (size_t i = 1; i < min(cf_coeffs.size(),
                             depth - offset_coeffs.size()); i++) {
    // initialize
    current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;
    
    prev_num2 = prev_num1;
    prev_num1 = current_num;
    
    prev_denom2= prev_denom1;
    prev_denom1 = current_denom;
    
    //rescale to avoid over- and underflow
    const double rescale_val = get_rescale_value(current_num, current_denom);
    
    current_num = current_num*rescale_val;
    current_denom = current_denom*rescale_val;
    
    prev_num1 = prev_num1*rescale_val;
    prev_num2 = prev_num2*rescale_val;

    prev_denom1 = prev_denom1*rescale_val;
    prev_denom2 = prev_denom2*rescale_val;
  }

  double offset_part = 0.0;
  for (size_t i = 0; i < offset_coeffs.size(); i++)
    offset_part += offset_coeffs[i]*pow(val, i);
  
  return offset_part + pow(val, min(depth, offset_coeffs.size()))*
    current_num/current_denom;
} 


// calculate ContinuedFraction approx when lower_offdiag > 0
static double 
evaluate_below_diagonal(const vector<double> &cf_coeffs,
                        const vector<double> &offset_coeffs,
                        const double val, const size_t depth) {
  
  //initialize
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;

  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0; 

  for (size_t i = 1; i < min(cf_coeffs.size(),
                             depth - offset_coeffs.size()); i++) {

    // recursion
    current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;

    prev_num2 = prev_num1;
    prev_num1 = current_num;

    prev_denom2= prev_denom1;
    prev_denom1 = current_denom;

    const double rescale_val = get_rescale_value(current_num, current_denom);

    current_num = current_num*rescale_val;
    current_denom = current_denom*rescale_val;

    prev_num1 = prev_num1*rescale_val;
    prev_num2 = prev_num2*rescale_val;
    
    prev_denom1 = prev_denom1*rescale_val;
    prev_denom2 = prev_denom2*rescale_val;
  }
  
  double offset_terms = 0.0;
  for (size_t i = 0; i < min(offset_coeffs.size(), depth); i++)
    offset_terms += offset_coeffs[i]*pow(val, i);
  
  // recall that if lower_offset > 0, we are working with 1/f, invert approx
  return 1.0/(offset_terms + pow(val, min(offset_coeffs.size(),depth))*
              current_num/current_denom);
}


// calculate ContinuedFraction approx when there is no offset
// uses euler's recursion
static double
evaluate_on_diagonal(const vector<double> &cf_coeffs, 
                     const double val, const size_t depth) {
  
  // initialize
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;

  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0; 

  for (size_t i = 1; i < min(cf_coeffs.size(), depth); i++) {
    // recursion
    current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;

    prev_num2 = prev_num1;
    prev_num1 = current_num;

    prev_denom2= prev_denom1;
    prev_denom1 = current_denom;

    const double rescale_val = get_rescale_value(current_num, current_denom);
    
    current_num = current_num*rescale_val;
    current_denom = current_denom*rescale_val;
    
    prev_num1 = prev_num1*rescale_val;
    prev_num2 = prev_num2*rescale_val;

    prev_denom1 = prev_denom1*rescale_val;
    prev_denom2 = prev_denom2*rescale_val;
  }
  return current_num/current_denom;
}


// calculate cont_frac approx depending on offset
double
ContinuedFraction::operator()(const double val) const {
  if (diagonal_idx > 0)
    return evaluate_above_diagonal(cf_coeffs, offset_coeffs, val, degree);
  
  if (diagonal_idx < 0)
    return evaluate_below_diagonal(cf_coeffs, offset_coeffs, val, degree);
  
  return evaluate_on_diagonal(cf_coeffs, val, degree);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
//////  COMPLEX NUMBER FUNCTIONS BELOW HERE
//////

// compute ContFrac_eval for complex values to compute deriv when no offset
static void
evaluate_complex_on_diagonal(const vector<double> &cf_coeffs,
                             const complex<double> perturbed_val,
                             const size_t depth, complex<double> &approx) {
  const complex<double> sqrt_neg1(0.0,1.0);
  if (norm(perturbed_val) == 0.0)
    approx = 0.0*sqrt_neg1;
  
  else {
    
    // Previous elements of the table to recursively fill it
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);
    
    for (size_t j = 1; j < min(cf_coeffs.size(), depth); j++) {
      //euler's recursion
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;

      //rescale to avoid over- and underflow
      const double rescale_val = get_rescale_value(current_num, current_denom);
      
      current_num = current_num*rescale_val;
      current_denom = current_denom*rescale_val;
      
      prev_num1 = prev_num1*rescale_val;
      prev_num2 = prev_num2*rescale_val;
      
      prev_denom1 = prev_denom1*rescale_val;
      prev_denom2 = prev_denom2*rescale_val;
    }
    approx = current_num/current_denom;
  }
}

// compute complex ContFrac_eval when above_diagonal > 0
static void
evaluate_complex_above_diagonal(const vector<double> &cf_coeffs,
                                const vector<double> &offset_coeffs,
                                const complex<double> perturbed_val,
                                const size_t depth, complex<double> &approx) {
  
  const complex<double> sqrt_neg1(0.0,1.0);
  if (norm(perturbed_val) == 0.0)
    approx = 0.0*sqrt_neg1;
  
  else {
    //initialize
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);

    for (size_t j = 1; j < min(depth - offset_coeffs.size(),
                               cf_coeffs.size()); j++) {
      
      //eulers recursion
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;

      //rescale to avoid over and underflow
      const double rescale_val = get_rescale_value(current_num, current_denom);
      
      current_num = current_num*rescale_val;
      current_denom = current_denom*rescale_val;
      
      prev_num1 = prev_num1*rescale_val;
      prev_num2 = prev_num2*rescale_val;
      
      prev_denom1 = prev_denom1*rescale_val;
      prev_denom2 = prev_denom2*rescale_val;
    }

    complex<double> offset_terms(0.0, 0.0);
    for (size_t i = 0; i < min(offset_coeffs.size(), depth); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);
    
    approx = 
      (offset_terms + pow(perturbed_val, min(offset_coeffs.size(), depth))*
       current_num/current_denom);
  }
} 

// compute cf approx when lower_offset > 0
static void
evaluate_complex_below_diagonal(const vector<double> &cf_coeffs,
                                const vector<double> &offset_coeffs,
                                const complex<double> perturbed_val,
                                const size_t depth,
                                complex<double> &approx) {
  const complex<double> sqrt_neg1(0.0,1.0);
  if (norm(perturbed_val) == 0.0)
    approx = 0.0*sqrt_neg1;
  else{
    // initialize
    complex<double> current_num(0.0, 0.0);
    complex<double> prev_num1(cf_coeffs[0], 0.0), prev_num2(0.0, 0.0);
    complex<double> current_denom(0.0, 0.0);
    complex<double> prev_denom1(1.0, 0.0), prev_denom2(1.0, 0.0);

    for (size_t j = 1; j < min(depth - offset_coeffs.size(),
                               cf_coeffs.size()); j++) {
      
      // euler's recursion
      complex<double> coeff(cf_coeffs[j], 0.0);
      current_num = prev_num1 + coeff*perturbed_val*prev_num2;
      current_denom = prev_denom1 + coeff*perturbed_val*prev_denom2;
      
      prev_num2 = prev_num1;
      prev_num1 = current_num;
      
      prev_denom2 = prev_denom1;
      prev_denom1 = current_denom;

      //rescale to avoid over and underflow
      const double rescale_val = get_rescale_value(current_num, current_denom);
      
      current_num = current_num*rescale_val;
      current_denom = current_denom*rescale_val;
      
      prev_num1 = prev_num1*rescale_val;
      prev_num2 = prev_num2*rescale_val;
      
      prev_denom1 = prev_denom1*rescale_val;
      prev_denom2 = prev_denom2*rescale_val;
    }
    
    complex<double> offset_terms(0.0, 0.0);
    for (size_t i = 0; i < min(offset_coeffs.size(), depth); i++)
      offset_terms += offset_coeffs[i]*pow(perturbed_val, i);

    approx = 1.0/
      (offset_terms + pow(perturbed_val, min(offset_coeffs.size(), depth))*
       current_num/current_denom);
  }
} 

/* compute cf approx for complex depending on offset df/dx =
 * lim_{delta -> 0} Imag(f(val+i*delta))/delta
 */
double
ContinuedFraction::complex_deriv(const double val) const {
  vector<double> ContFracCoeffs(cf_coeffs);
  vector<double> ContFracOffCoeffs(offset_coeffs);
  
  const complex<double> sqrt_neg1(0.0,1.0);
  complex<double> df(0.0, 0.0);
  complex<double> value(val, 0.0);
  
  if (diagonal_idx == 0)
    evaluate_complex_on_diagonal(ContFracCoeffs, value + DERIV_DELTA*sqrt_neg1, degree, df);
  
  else if (diagonal_idx > 0)
    evaluate_complex_above_diagonal(ContFracCoeffs, ContFracOffCoeffs,
                                    value + DERIV_DELTA*sqrt_neg1, degree, df);
  
  else if (diagonal_idx < 0)
    evaluate_complex_below_diagonal(ContFracCoeffs, ContFracOffCoeffs,
                                    value + DERIV_DELTA*sqrt_neg1, degree, df);
  
  return imag(df)/DERIV_DELTA;
}


std::ostream&
operator<<(std::ostream& the_stream, const ContinuedFraction &cf) {
  std::ios_base::fmtflags orig_flags = the_stream.flags();
  the_stream.setf(std::ios_base::fixed, std::ios_base::floatfield);
  the_stream.precision(2);
  the_stream << "OFFSET_COEFFS" << '\t' << "PS_COEFFS" << '\n';
  const size_t offset = cf.offset_coeffs.size();
  for (size_t i = 0; i < offset; ++i)
    the_stream << std::setw(12) << cf.offset_coeffs[i] << '\t'
               << std::setw(12) << cf.ps_coeffs[i] << '\n';
  the_stream << "CF_COEFFS" << '\n';
  for (size_t i = 0; i < cf.cf_coeffs.size(); ++i)
    the_stream << std::setw(12) << cf.cf_coeffs[i] << '\t'
               << std::setw(12) << cf.ps_coeffs[i + offset] << '\n';
  the_stream.flags(orig_flags);
  return the_stream;
}

// Extrapolates the curve, for given values (step & max) and numbers
// of terms
void
ContinuedFraction::extrapolate_distinct(const vector<double> &counts_hist,
                                        const double max_value, 
                                        const double step_size,
                                        vector<double> &estimates) const {
  const double hist_sum = accumulate(counts_hist.begin(), counts_hist.end(), 0.0);
  estimates.clear();
  estimates.push_back(hist_sum);
  for (double t = step_size; t <= max_value; t += step_size)
    estimates.push_back(hist_sum + t*operator()(t));
}


void
ContinuedFraction::extrapolate_count(const vector<double> &counts_hist,
                                     const double max_value,
                                     const double step_size,
                                     const size_t count,
                                     vector<double> &estimates) const {
  const double current_count = counts_hist[count];
  estimates.clear();
  estimates.push_back(current_count);
  for (double t = step_size; t <= max_value; t += step_size)
    estimates.push_back(operator()(t));
}


void
ContinuedFraction::extrapolate_mincount(const vector<double> &counts_hist,
                                        const double max_value,
                                        const double step_size,
                                        const size_t mincount,
                                        vector<double> &estimates) const {
  const double current_observed = 
    accumulate(counts_hist.begin() + mincount, counts_hist.end(), 0.0);
  //  cerr << "current_observed = " << current_observed << endl;

  estimates.clear();
  estimates.push_back(current_observed);
  for (double t = step_size; t <= max_value; t += step_size)
    estimates.push_back(current_observed + t*operator()(t));
}

void
ContinuedFraction::extrapolate_saturation(const vector<double> &counts_hist,
                                          const double vals_sum,
					  const double initial_extrap_val,
                                          const double max_value,
                                          const double step_size,
                                          vector<double> &saturation_estimates) const {
  for(double t = initial_extrap_val; t < max_value; t += step_size)
    saturation_estimates.push_back(operator()(t)/vals_sum);
}

void
ContinuedFraction::extrapolate_yield_deriv(const vector<double> &counts_hist,
                                           const double vals_sum,
					   const double initial_extrap_val, 
                                           const double max_value,
                                           const double step_size,
                                           vector<double> &saturation_estimates) const {
  //  saturation_estimates.push_back(counts_hist[1]/vals_sum);
  for(double t = initial_extrap_val; t < max_value; t += step_size)
    saturation_estimates.push_back((operator()(t) 
                                    + t*complex_deriv(t))/vals_sum);
}


//////////////////////////////////////////////////////////
// Y50: expected # reads to have 50% distinct (or 50% duplicates)
// A measure of library quality 
static double
sample_count_distinct(const gsl_rng *rng,
		      const vector<size_t> &full_umis,
		      const size_t sample_size) {
  vector<size_t> sample_umis(sample_size);
  gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
		 (size_t *)&full_umis.front(), full_umis.size(), 
		 sizeof(size_t));
  double count = 1.0;
  for (size_t i = 1; i < sample_umis.size(); i++)
    if(sample_umis[i] != sample_umis[i-1])
      count++;

  return count;
}


// calculate the expected number of reads to reach
// dupl_level% saturation
// use bisection since the yield curve is concave
// assuming the CF is optimal
double
ContinuedFraction::Ylevel(const vector<double> &counts_hist, const double dupl_level,
			  const double vals_sum, const double max_val,
			  const double tol, const size_t max_iter) const {

  const double observed_distinct = 
    accumulate(counts_hist.begin(), counts_hist.end(), 0.0);

  // case 1: the observed library is already above dupl_level% duplicates
    // search by bisection, subsampling the library
  if(observed_distinct < dupl_level*vals_sum){
    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());

    // set to sample from
    vector<size_t> full_umis;
    size_t indx = 1;
    for (size_t i = 1; i < counts_hist.size(); i++){
      for (size_t j = 0; j < counts_hist[i]; j++){
	for (size_t k = 0; k < i; k++){
	  full_umis.push_back(indx);
	}
	indx++;
      }
    }

    double upper_distinct = observed_distinct;
    size_t upper_sample_size = static_cast<size_t>(vals_sum);
    double lower_distinct = 0.0;
    size_t lower_sample_size = 0;
    size_t mid_sample_size = upper_sample_size/2;
    double mid_distinct = sample_count_distinct(rng, full_umis, mid_sample_size);
    
    size_t iter = 0;
    double rel_error = fabs(mid_distinct - dupl_level*mid_sample_size)/mid_sample_size;
    while(rel_error > tol && iter < max_iter){
      // if observed_distinct < dupl_level*sample_size, the intersection is lower
      if(mid_distinct < dupl_level*mid_sample_size){
	upper_sample_size = mid_sample_size;
	upper_distinct = mid_distinct;
	mid_sample_size = (upper_sample_size + lower_sample_size)/2;
	mid_distinct = sample_count_distinct(rng, full_umis, mid_sample_size);
      }
      // if observed_distinct > dupl_level*sample_size, the intersection is higher
      else if(mid_distinct > dupl_level*mid_sample_size){
	lower_sample_size = mid_sample_size;
	lower_distinct = mid_distinct;
	mid_sample_size = (upper_sample_size + lower_sample_size)/2;
	mid_distinct = sample_count_distinct(rng, full_umis, mid_sample_size);
      }

      rel_error = fabs(mid_distinct - dupl_level*mid_sample_size)/mid_sample_size;
      iter++;
    }
    // return estimated sample size as double
    return static_cast<double>(mid_sample_size);
  }

  // case 2: observed distinct is less that dupl_level% of sample size
  // need to extrapolate
  else{
    double upper_val = max_val;
    double upper_distinct = observed_distinct + upper_val*operator()(upper_val);
    double upper_sample_size = vals_sum*(upper_val + 1.0);
    double lower_val = 0.0;
    double lower_distinct = observed_distinct;
    double lower_sample_size = vals_sum;

    // max_val to low, double it
    while(upper_distinct > dupl_level*upper_sample_size){
      lower_val = upper_val;
      lower_distinct = upper_distinct;
      lower_sample_size = upper_sample_size;
      upper_val = 2.0*upper_val;
      upper_distinct = observed_distinct + upper_val*operator()(upper_val);
      upper_sample_size = vals_sum*(upper_val + 1.0);
    }

    double mid_val = (upper_val + lower_val)/2.0;
    double mid_distinct = observed_distinct + mid_val*operator()(mid_val);
    double mid_sample_size = vals_sum*(mid_val + 1.0);

    // find Y50 by bisection
    size_t iter = 0;
    double rel_error = fabs(mid_distinct - dupl_level*mid_sample_size)/mid_sample_size;
    while(rel_error > tol && iter < max_iter){
      // if observed_distinct < dupl_level*sample_size, the intersection is lower
      if(mid_distinct < dupl_level*mid_sample_size){
	upper_val = mid_val;
	upper_sample_size = mid_sample_size;
	upper_distinct = mid_distinct;
	mid_val = (upper_val + lower_val)/2.0;
	mid_sample_size = vals_sum*(mid_val + 1.0);
	mid_distinct = observed_distinct + mid_val*operator()(mid_val);
      }
      // if observed_distinct > dupl_level*sample_size, the intersection is higher
      else if(mid_distinct > dupl_level*mid_sample_size){
	lower_val = mid_val;
	lower_sample_size = mid_sample_size;
	lower_distinct = mid_distinct;
	mid_val = (upper_val + lower_val)/2.0;
	mid_sample_size = vals_sum*(mid_val + 1.0);
	mid_distinct = observed_distinct + mid_val*operator()(mid_val);
      }

      rel_error = fabs(mid_distinct - dupl_level*mid_sample_size)/mid_sample_size;
      iter++;
    }

    return mid_sample_size;
  }

}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////
/////////////////
////////////////  CONTINUED FRACTION APPROXIMATION CLASS BELOW
///////////////
//////////////
/////////////
////////////


// calculate cf_coeffs depending on offset
ContinuedFractionApproximation::ContinuedFractionApproximation(const int di, const size_t mt, 
                                                               const double ss, const double mv) :
  diagonal_idx(di), max_terms(mt), step_size(ss), max_value(mv) {}


static inline double
movement(const double a, const double b) {
  return fabs((a - b)/std::max(a, b)); //delta
}

/* locate zero deriv by bisection to find local max within (prev_val,
   val)
*/ 
// need to modify to account for change in evaluating
double
ContinuedFractionApproximation::locate_zero_cf_deriv(const ContinuedFraction &cf, 
                                                     const double val, 
                                                     const double prev_val) const {
  
  double val_low = prev_val;
  double deriv_low = cf.complex_deriv(val_low);
  
  double val_high = val;
  double deriv_high = cf.complex_deriv(val_high);
  
  double val_mid = (val - prev_val)/2.0;
  double deriv_mid = std::numeric_limits<double>::max();
  
  double diff = std::numeric_limits<double>::max();
  double prev_deriv = std::numeric_limits<double>::max();

  while (diff > TOLERANCE && movement(val_low, val_high) > TOLERANCE) {
    
    val_mid = (val_low + val_high)/2.0;
    deriv_mid = cf.complex_deriv(val_mid);
    
    if ((deriv_mid > 0 && deriv_low < 0) || (deriv_mid < 0 && deriv_low > 0)) {
      val_high = val_mid;
      deriv_high = deriv_mid;
    }
    else {
      val_low = val_mid;
      deriv_low = deriv_mid;
    }
    diff = fabs((prev_deriv - deriv_mid)/prev_deriv);
    prev_deriv = deriv_mid;
  }
  
  return val_mid;
}


// search (min_val, max_val) for local max
// return location of local max
double
ContinuedFractionApproximation::local_max(const ContinuedFraction &cf,
                                          const double deriv_upper) const {
  double current_max = cf(0.0);
  for (double val = 100*SEARCH_STEP_SIZE; val <= SEARCH_MAX_VAL; val += 100*SEARCH_STEP_SIZE)
    current_max = std::max(current_max, cf(locate_zero_cf_deriv(cf, val, val - 100*SEARCH_STEP_SIZE)));
  return current_max;
}


/* Checks if estimates are stable (derivative large) for the
 * particular approximation (degrees of num and denom) at a specific
 * point
 */
static bool
check_yield_estimates_stability(const vector<double> &estimates) {
  // make sure that the estimate is increasing in the time_step and
  // is below the initial distinct per step_size
  for (size_t i = 1; i < estimates.size(); ++i){
    if (estimates[i] < estimates[i - 1] ){
      return false;
    }
    if(i >= 2 && (estimates[i] - estimates[i - 1] >
                  estimates[i - 1] - estimates[i - 2])){
      return false;
    }
  }
    

  // fake check
  /*  for(size_t i = 1; i < estimates.size(); ++i)
      if(estimates[i] < 0.0 || estimates[i] > 1e9)
      return false;
  */
  return true;
}


/* Finds the optimal number of terms (i.e. degree, depth, etc.) of the
 * continued fraction by checking for stability of estimates at
 * specific points for yield.
 */
// New way for searching for optimal CF
ContinuedFraction
ContinuedFractionApproximation::optimal_cont_frac_distinct(const vector<double> &counts_hist) const {
  //do this outside
  // ensure that we will use an underestimate
  //  const size_t local_max_terms = max_terms - (max_terms % 2 == 1); 
 
  if(max_terms >= counts_hist.size()){
    cerr << "max terms = " << max_terms << endl;
    cerr << "hist size = " << counts_hist.size() << endl;
  } 
  assert(max_terms < counts_hist.size());
  
  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    counts_sum += i*counts_hist[i];
  
  vector<double> full_ps_coeffs;
  for (size_t j = 1; j <= max_terms; j++)
    full_ps_coeffs.push_back(counts_hist[j]*pow(-1, j + 1));

  ContinuedFraction full_CF(full_ps_coeffs, -1, max_terms);  

  // if max terms = 4, check only that degree
  if(max_terms == 4 || max_terms == 3 
     || max_terms == 5 || max_terms == 6){   
    vector<double> estimates;
    full_CF.extrapolate_distinct(counts_hist, SEARCH_MAX_VAL, SEARCH_STEP_SIZE, estimates);
    // return the continued fraction if it is stable
    if (check_yield_estimates_stability(estimates))
      return full_CF;
  }
  else{
    //if max terms >= 8, start at 8 and check increasing cont frac's
    size_t curr_terms = 0;
    if(max_terms % 2 == 0)
      curr_terms = 8;
    else
      curr_terms = 7;
    while (curr_terms <= max_terms) {    
      ContinuedFraction curr_cf 
	= ContinuedFraction::truncate_degree(full_CF, curr_terms);
      vector<double> estimates;
      curr_cf.extrapolate_distinct(counts_hist, SEARCH_MAX_VAL, SEARCH_STEP_SIZE, estimates);
          
    // return the continued fraction if it is stable
      if (check_yield_estimates_stability(estimates))
	return curr_cf;
    
      curr_terms += 2;
    // if not cf not acceptable, increase degree
    }
  }
   // no stable continued fraction: return null
  return ContinuedFraction();  
}


/* Finds the optimal number of terms (i.e. degree, depth, etc.) of the
 * continued fraction by checking for stability of estimates at
 * specific points for yield.
 */
/*  The old way of searching for the optimal CF
ContinuedFraction
ContinuedFractionApproximation::optimal_cont_frac_distinct(const vector<double> &counts_hist) const {
  //do this outside
  // ensure that we will use an underestimate
  //  const size_t local_max_terms = max_terms - (max_terms % 2 == 1); 
  
  assert(max_terms < counts_hist.size());
  
  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    counts_sum += i*counts_hist[i];
  
  vector<double> ps_coeffs;

  for (size_t j = 1; j <= max_terms; j++)
    ps_coeffs.push_back(counts_hist[j]*pow(-1, j + 1));  

  ContinuedFraction curr_cf(ps_coeffs, -1, max_terms - 1);
  
  while (curr_cf.degree >= MIN_ALLOWED_DEGREE) {    
    // compute the estimates for the desired set of points
    vector<double> estimates;
    curr_cf.extrapolate_distinct(counts_hist, SEARCH_MAX_VAL, SEARCH_STEP_SIZE, estimates);
    
    // return the continued fraction if it is stable
    if (check_yield_estimates_stability(estimates))
      return curr_cf;
    
    // if not cf not acceptable, decrease degree
    curr_cf = ContinuedFraction::decrease_degree(curr_cf, 2);
  }
  
  //  throw SMITHLABException("unable to fit continued fraction");
  
  // no stable continued fraction: return null
  return ContinuedFraction();  
}
*/


/* Checks if count estimates do not go out of bounds 
 * and are unimodal
 */
static bool
check_count_estimates_stability(const vector<double> &estimates,
                                const size_t count) {
  size_t number_modes = 0;
  for(size_t i = 1; i < estimates.size(); ++i){
    // make sure estimates are in bounds
    if(estimates[i] < 0.0 || estimates[i] > 3.2e9/count || !finite(estimates[i]))
      return false;
    // count modes by detecting change in sign of derivative
    if(i < estimates.size() - 1 &&
       ((estimates[i] - estimates[i - 1])
        *(estimates[i + 1] - estimates[i]) < 0))
      number_modes++;
  }
  // check unimodality
  if(number_modes > 1)
    return false;
 
  return true;
}


static void
construct_count_ps_coeffs(const vector<double> &counts_hist,
                          const size_t max_terms, const size_t count,
                          vector<double> &ps_coeffs) {
  ps_coeffs.clear();
  for(size_t i = 0; i < max_terms; i++){
    double curr_coeff = 0.0;
    for(size_t l = 0; l <= i; l++){
      const double first_binom_coeff =
        exp(gsl_sf_lnfact(l + count) - gsl_sf_lnfact(l) - gsl_sf_lnfact(count));
      const double second_binom_coeff = 
        exp(gsl_sf_lnfact(count) - gsl_sf_lnfact(i - l) - gsl_sf_lnfact(count - i + l));
      curr_coeff += 
        pow(-1.0, l)*first_binom_coeff*counts_hist[count + l]*second_binom_coeff;
    }
    ps_coeffs.push_back(curr_coeff);
  }
}

ContinuedFraction
ContinuedFractionApproximation::optimal_cont_frac_count(const vector<double> &counts_hist,
                                                        const size_t count) const {
  

  
  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    counts_sum += i*counts_hist[i];
  
  vector<double> ps_coeffs;
  construct_count_ps_coeffs(counts_hist, max_terms, count, ps_coeffs);

  
 // for (size_t j = 0; j < max_terms; j++){
 //   const double binom_coeff =
 //     exp(gsl_sf_lnfact(j + count) - gsl_sf_lnfact(j) - gsl_sf_lnfact(count));
 //   ps_coeffs.push_back(pow(-1.0, j)*binom_coeff*counts_hist[j + count]);  
 // }
  int order = 0;

  if(max_terms < MIN_ALLOWED_DEGREE + abs(order)){
    cerr << "max_terms too small" << endl;
    return ContinuedFraction();
  }

  ContinuedFraction curr_cf(ps_coeffs, order, max_terms);
  
  while (curr_cf.degree >= MIN_ALLOWED_DEGREE + abs(order)) {    
    // compute the estimates for the desired set of points
    vector<double> estimates;
    curr_cf.extrapolate_count(counts_hist, SEARCH_MAX_VAL, SEARCH_STEP_SIZE, count, estimates);
    
    // return the continued fraction if it is stable
    if (check_count_estimates_stability(estimates, count))
      return curr_cf;

    
    // if not cf not acceptable, decrease degree
    curr_cf = ContinuedFraction::decrease_degree(curr_cf, 2);
  }
  
  //  throw SMITHLABException("unable to fit continued fraction");
  
  // no stable continued fraction: return null
  return ContinuedFraction();  
}
  

/*
static void
construct_mincount_ps_coeffs(const vector<double> &counts_hist,
                             const size_t max_terms, const size_t mincount,
                             vector<double> &ps_coeffs) {
  ps_coeffs.clear();
  for(size_t i = 1; i < max_terms; i++){
    double curr_coeff = 0.0;
    for(size_t r = 1; r < mincount; r++){
      for(size_t l = 0; l <= i; l++){
        const double first_binom_coeff =
          exp(gsl_sf_lnfact(l + r) - gsl_sf_lnfact(l) - gsl_sf_lnfact(r));
        const double second_binom_coeff = 
          exp(gsl_sf_lnfact(r) - gsl_sf_lnfact(i - l) - gsl_sf_lnfact(r - i + l));
        curr_coeff += 
          pow(-1.0, l + 1)*first_binom_coeff*counts_hist[r + l]*second_binom_coeff;
      }
    }
    ps_coeffs.push_back(curr_coeff);
  }
}

static bool
check_mincount_estimates_stability(const vector<double> &estimates,
                                   const double max_change_per_time_step) {
  // make sure that the estimate is increasing in the time_step and
  // is below the initial distinct per step_size
  for (size_t i = 1; i < estimates.size(); ++i){
    if(!finite(estimates[i])){
      return false;
    }
    if ((estimates[i] < estimates[i - 1]) ||
        (estimates[i] - estimates[i - 1] > max_change_per_time_step)){
      return false;
    }
  }
    

  return true;
}


ContinuedFraction
ContinuedFractionApproximation::optimal_cont_frac_mincount(const vector<double> &counts_hist,
                                                           const size_t mincount) const {
  

  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    counts_sum += i*counts_hist[i];

  double observed_mincount = accumulate(counts_hist.begin() + mincount, counts_hist.end(), 0.0);
  
  vector<double> ps_coeffs;
  construct_mincount_ps_coeffs(counts_hist, max_terms, mincount, ps_coeffs);

  const int distinct_order = -1;
  const int remainder_order = 0;

  // if max_terms is too small, unacceptable extrapolation
  if(max_terms < MIN_ALLOWED_DEGREE + abs(order))
    return ContinuedFraction();
  

  ContinuedFraction curr_cf(ps_coeffs, order, max_terms);

  while (curr_cf.degree >= MIN_ALLOWED_DEGREE + abs(order)) {    
    // compute the estimates for the desired set of points
    vector<double> estimates;

    curr_cf.extrapolate_mincount(counts_hist, SEARCH_MAX_VAL, 
                                 SEARCH_STEP_SIZE, mincount, estimates);
    const double max_reads_per_step = counts_sum*SEARCH_STEP_SIZE;
    // return the continued fraction if it is stable
    if (check_mincount_estimates_stability(estimates, max_reads_per_step))
      return curr_cf;

    
    // if not cf not acceptable, decrease degree
    curr_cf = ContinuedFraction::decrease_degree(curr_cf, 2);
  }
  
  //  throw SMITHLABException("unable to fit continued fraction");
  
  // no stable continued fraction: return null
  return ContinuedFraction();  
}
*/



static void
construct_mincount_ps_coeffs(const vector<double> &counts_hist,
                             const size_t max_terms, const size_t mincount,
                             vector<double> &ps_coeffs) {
  ps_coeffs.clear();
  for(size_t i = 1; i <= max_terms; i++){
    double curr_coeff = 0.0;
    for(size_t l = 0; l < mincount; l++){
      for(size_t j = 0; j <= i; j++){
	if((l >= i - j) && (l + j > 0)){
          const double binom_coeff_1 = 
	    exp(gsl_sf_lnfact(j + l) - gsl_sf_lnfact(j)
		-gsl_sf_lnfact(l));
	  const double binom_coeff_2 = 
	    exp(gsl_sf_lnfact(l) - gsl_sf_lnfact(i - j)
		-gsl_sf_lnfact(l - (i - j)));
	  curr_coeff += 
	    pow(-1, j + 1)*binom_coeff_1*binom_coeff_2*counts_hist[l + j];
	}
      }
    }
    ps_coeffs.push_back(curr_coeff);
  }
}

static bool
check_mincount_estimates_stability(const vector<double> &estimates,
                                   const double max_change_per_time_step) {
  // make sure that the estimate is increasing in the time_step and
  // is below the initial distinct per step_size
  for (size_t i = 1; i < estimates.size(); ++i){
    if(!finite(estimates[i])){
      //    cerr << "not finite at " << i << "\t" << estimates[i] << endl;
      return false;
    }
    if ((estimates[i] < estimates[i - 1]) ||
        (estimates[i] - estimates[i - 1] > max_change_per_time_step)){
      /*     if(estimates[i] < estimates[i - 1])
	cerr << i - 1 << "th estimate (" << estimates[i-1] 
	     << ") greater than " << i << "th (" << estimates[i] << ")" << endl;
      if(estimates[i] - estimates[i - 1] > max_change_per_time_step)
	cerr << "change is more than max_change_per_time_step, change = "
	     << estimates[i] - estimates[i - 1] 
	     << ", max_change_per_time_step = " << max_change_per_time_step << endl;
      */
      return false;
    }
  }
    
  return true;
}


ContinuedFraction
ContinuedFractionApproximation::optimal_cont_frac_mincount(const vector<double> &counts_hist,
                                                           const size_t mincount,
const int order) const {
  

  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    counts_sum += i*counts_hist[i];
  const double max_reads_per_step = counts_sum*SEARCH_STEP_SIZE;

  vector<double> ps_coeffs;
  construct_mincount_ps_coeffs(counts_hist, max_terms, mincount, ps_coeffs);

  /*
  cerr << "ps_coeffs = " << endl;
  for(size_t i = 0; i < ps_coeffs.size(); i++)
    cerr << ps_coeffs[i] << endl;
  cerr << endl;
  */

  // if max_terms is too small, unacceptable extrapolation
  if(max_terms < MIN_ALLOWED_DEGREE + abs(order))
    return ContinuedFraction();
  

  ContinuedFraction curr_cf(ps_coeffs, order, max_terms);

  while (curr_cf.degree >= MIN_ALLOWED_DEGREE + abs(order)) {    
    // compute the estimates for the desired set of points
    vector<double> estimates;

    curr_cf.extrapolate_mincount(counts_hist, SEARCH_MAX_VAL, 
                                 SEARCH_STEP_SIZE, mincount, estimates);
    // return the continued fraction if it is stable
    if (check_mincount_estimates_stability(estimates, max_reads_per_step))
      return curr_cf;

    
    // if not cf not acceptable, decrease degree
    curr_cf = ContinuedFraction::decrease_degree(curr_cf, 2);
  }
  
  //  throw SMITHLABException("unable to fit continued fraction");
  
  // no stable continued fraction: return null
  return ContinuedFraction();  
}

static inline bool
check_saturation_estimates(const vector<double> estimates){
  if(estimates.empty())
    return false;

  // make sure estimates are decreasing and
  // between 0 & 1
  if(estimates[0] >= 1.0 || estimates[0] < 0.0)
    return false;

  for(size_t i = 1; i < estimates.size(); i++)
    if(estimates[i] > estimates[i-1] ||
       estimates[i] >= 1.0 ||
       estimates[i] < 0.0) 
      return false;
  
  return true;
}

/* Finds the optimal number of terms (i.e. degree, depth, etc.) of the
 * continued fraction by checking for stability of estimates at
 * specific points for saturation estimates.
 */
ContinuedFraction
ContinuedFractionApproximation::optimal_cont_frac_satur(const vector<double> &counts_hist) const { 

  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    counts_sum += i*counts_hist[i];
  
  vector<double> full_ps_coeffs;
  for (size_t j = 1; j < max_terms; j++)
    full_ps_coeffs.push_back(counts_hist[j]*j*pow(-1, j + 1));

  ContinuedFraction full_CF(full_ps_coeffs, -1, max_terms); 

  if(max_terms == 4){
    vector<double> estimates;
    full_CF.extrapolate_yield_deriv(counts_hist, counts_sum, 0.0,  
				    SEARCH_MAX_VAL, SEARCH_STEP_SIZE, 
				    estimates);
    if (check_saturation_estimates(estimates))
      return full_CF;
  }
  else{
    size_t curr_terms = 6;
    while (curr_terms <= max_terms){
     ContinuedFraction curr_cf 
	= ContinuedFraction::truncate_degree(full_CF, curr_terms);
      vector<double> estimates;
      curr_cf.extrapolate_saturation(counts_hist, counts_sum, 0.0,  
				      SEARCH_MAX_VAL, SEARCH_STEP_SIZE, 
				      estimates);
    if (check_saturation_estimates(estimates))
      return curr_cf;

    curr_terms += 2;
    }
  }

  // no stable continued fraction: return null
  return ContinuedFraction();  
}
