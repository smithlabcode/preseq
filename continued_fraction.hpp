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

#ifndef CONTINUED_FRACTION_HPP
#define CONTINUED_FRACTION_HPP

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <complex>
#include <cassert>
#include <unistd.h>



struct ContinuedFraction {
  // Constructors
  ContinuedFraction() {}
  ContinuedFraction(const std::vector<double> &ps_cf, 
                    const int di, const size_t dg);

  // Evaluate the continued fraction
  double operator()(const double val) const;

  //////////////////////////////////////////
  // Extrapolation functions

  // Evaluate the continued fraction estimating distinct
  // along a curve from 0 to max_value
  void 
  extrapolate_distinct(const std::vector<double> &counts_hist,
                       const double max_value, const double step_size,
                       std::vector<double> &estimates) const;
  
  
  bool is_valid() const {return !cf_coeffs.empty();}
  size_t return_degree() const {return degree;}


  static ContinuedFraction truncate_degree(const ContinuedFraction &fullCF,
					   const size_t truncated_degree);
  
  std::vector<double> ps_coeffs;
  std::vector<double> cf_coeffs;
  std::vector<double> offset_coeffs;
  int diagonal_idx;
  size_t degree;
};

std::ostream& 
operator<<(std::ostream& the_stream, const ContinuedFraction &cf);


class ContinuedFractionApproximation {
public:
  // Constructor
  ContinuedFractionApproximation(const int di, const size_t mt, 
                                 const double ss, const double mv);
  
  //find best cont frac approx for estimating distinct
  ContinuedFraction
  optimal_cont_frac_distinct(const std::vector<double> &counts_hist) const;

private:
  
  int diagonal_idx; // the diagonal to work with for estimates
  size_t max_terms; // the maximum number of terms to try for a CF
  double step_size; // the step size to use when training
  double max_value; // the largest value to check when training

  static const size_t MIN_ALLOWED_DEGREE = 4;
  
  // largest value to search for lowerbound and stability
  static const double SEARCH_MAX_VAL = 200; 
  
  //step size for search of lowerbound and stability
  static const double SEARCH_STEP_SIZE = 0.01; 

};

#endif
 
