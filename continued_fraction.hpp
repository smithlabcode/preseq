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

#include <numeric>
#include <vector>
#include <fstream>
#include <iomanip>



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
  extrapolate_distinct(const double initial_sum,
                       const double max_value, const double step_size,
                       std::vector<double> &estimates) const;

  void
  extrapolate_distinct(const double max_value, const double step_size,
                       std::vector<double> &estimates) const;
  
  bool is_valid() const {return !cf_coeffs.empty();}
  size_t return_degree() const {return degree;}

  // Return new ContinuedFraction with degree decrement less than CF
  static ContinuedFraction decrease_degree(const ContinuedFraction &CF,
                                           const size_t decrement);

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
  ContinuedFractionApproximation(const int di, const size_t mt);
  
  //find best cont frac approx for estimating distinct
  ContinuedFraction
  optimal_cont_frac_distinct(const std::vector<double> &counts_hist) const;

  int get_diagonal() const {return diagonal_idx;}

  //find best cont frac approx for a power series
  ContinuedFraction
  optimal_powerseries_to_cont_frac(const std::vector<double> &ps_coeff,
                                   const double step_size) const;

private:
  
  int diagonal_idx; // the diagonal to work with for estimates
  size_t max_terms; // the maximum number of terms to try for a CF

  static const size_t MIN_ALLOWED_DEGREE;
  
  // largest value to search for lowerbound and stability
  static const double SEARCH_MAX_VAL; 
  
  //step size for search of lowerbound and stability
  static const double SEARCH_STEP_SIZE; 

};

#endif
