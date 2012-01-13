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

struct cont_frac {
  // Constructor
  cont_frac(const std::vector<double> &coeffs,
		    const size_t in_lower, const size_t in_upper) :
    lower_offset(in_lower), upper_offset(in_upper), ps_coeffs(coeffs)
  {;}

  // Mutators
  void set_lower_offset(const size_t l_o) {lower_offset = l_o; 
    upper_offset = 0;}
  void set_upper_offset(const size_t u_o) {upper_offset = u_o; 
    lower_offset = 0;}
  void set_ps_coeffs(const std::vector<double> &coeffs)
  {ps_coeffs = coeffs;}
  void set_cf_coeffs(const std::vector<double> &coeffs)
  {cf_coeffs = coeffs;}
  void set_offset_coeffs(const std::vector<double> &coeffs)
  {offset_coeffs = coeffs;}
  
  // Accessors
  void get_cf_coeffs(std::vector<double> &return_coeffs) const
  {return_coeffs = cf_coeffs;}
  void get_offset_coeffs(std::vector<double> &return_coeffs) const 
  {return_coeffs = offset_coeffs;}
  void get_ps_coeffs(std::vector<double> &return_coeffs) const
  {return_coeffs = ps_coeffs;}
  
  // Evaluators
  double evaluate(const double val, const size_t depth);
  double complex_deriv(const double val, const size_t depth);

  size_t lower_offset;
  size_t upper_offset;
  std::vector<double> ps_coeffs;
  std::vector<double> cf_coeffs;
  std::vector<double> offset_coeffs;
};

class ContFracApprox {
public:
  static const size_t MINIMUM_ALLOWED_DEGREE = 6;

  // Constructor
  ContFracApprox(const cont_frac &cf_instance, const size_t max_terms) :
    cont_frac_estimate(cf_instance), depth(max_terms) {compute_cf_coeffs();}

  // Mutators
  void compute_cf_coeffs();
  void set_ps_coeffs(std::vector<double> &coeffs)
  {cont_frac_estimate.set_ps_coeffs(coeffs); compute_cf_coeffs();}
  void set_depth(const size_t max_terms); 

  // estimators
  double evaluate(const double val) 
  {return cont_frac_estimate.evaluate(val, depth);}
  double complex_deriv(const double val)
  {return cont_frac_estimate.evaluate(val, depth);}
  double locate_local_max(const double lower_limit, 
			  const double upper_limit,
			  const double step_size, 
			  const double upper_bound,
			  const double deriv_upper_bound);


private:
  cont_frac cont_frac_estimate;
  size_t depth;

  double locate_zero_cf_deriv(const double val, const double prev_val);
};


#endif
