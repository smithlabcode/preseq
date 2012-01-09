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

#ifndef CONTINUED_FRACTION_HPP
#define CONTINUED_FRACTION_HPP

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <complex>

class ContinuedFraction {
public:
  // Constructor
  ContinuedFraction(const std::vector<double> in_cfs, 
		    const std::vector<double> in_offset,
		    const size_t in_lower, const size_t in_upper) :
    cf_coeffs(in_cfs), offset_coeffs(in_offset), lower_offset(in_lower), 
    upper_offset(in_upper) {}

  // Mutators
  void set_lower_offset(const size_t l_o) {lower_offset = l_o; upper_offset = 0;}
  void set_upper_offset(const size_t u_o) {upper_offset = u_o; lower_offset = 0;}
  void set_offset_coeffs(const std::vector<double> coeffs) {offset_coeffs = coeffs;}
  void set_cf_coeffs(const std::vector<double> coeffs) {cf_coeffs = coeffs;}
  
  // Accessors
  void get_cf_coeffs(std::vector<double> &return_coeffs) const
  {return_coeffs = cf_coeffs;}
  void get_offset_coeffs(std::vector<double> &return_coeffs) const 
  {return_coeffs = offset_coeffs;}
  
  // Push this inside
  void compute_cf_coeffs(const std::vector<double> ps_coeffs, const size_t depth) const;
  double cf_approx(const double time, const double tolerance);
  double cf_deriv_complex(const double val, const double dx,
                          const double tolerance);
  
  // Should not be public -- functionality should be pushed inside
  double locate_zero_cf_deriv(const double val, const double prev_val,
                              const double dx, const double tolerance);
  
private:
  std::vector<double> cf_coeffs;
  std::vector<double> offset_coeffs;
  size_t lower_offset;
  size_t upper_offset;
};

// double
// evaluate_continued_fraction(const std::vector<double> &cf_coeffs, const double t);

#endif
