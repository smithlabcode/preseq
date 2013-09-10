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
  vector< vector<double> > q_table(depth, vector<double>(depth+1, 0.0)); //an array that is depth x depth+1
  vector< vector<double> > e_table(depth, vector<double>(depth+1, 0.0));

  for (size_t j = 0; j < depth-1; j++) //fill first position in the q_table with a vector that holds the ratio of power series coeff Cn to Cn-1
    q_table[1][j] = ps_coeffs[j + 1]/ps_coeffs[j];
  
  for (size_t j = 0; j < depth-1; j++) //fill first position in e_table with the quotient difference algorithm relation 
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
      //alternate appending the value in e[i][0] and value in q[i][0] (change signs) to the vector for cont. frac. coeffs
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

ContinuedFraction::ContinuedFraction(const vector<double> &ps_cf, 
                                     const int di, const size_t dg) :
  ps_coeffs(ps_cf), diagonal_idx(di), degree(dg) {

  if (diagonal_idx == 0)
    quotdiff_algorithm(ps_coeffs, cf_coeffs);
  else if (diagonal_idx > 0)
    quotdiff_above_diagonal(ps_coeffs, diagonal_idx, cf_coeffs, offset_coeffs);
  else // if(cont_frac_estimate.lower_offset > 0) {
    quotdiff_below_diagonal(ps_coeffs, -diagonal_idx, cf_coeffs, offset_coeffs);
  // notice the "-" above...
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
  for (size_t i = 0; i < offset_coeffs.size(); i++)
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

  for (size_t j = 1; j < max_terms; j++)
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
