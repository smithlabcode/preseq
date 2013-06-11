/*    ZTNB:
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Timothy Daley
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


#ifndef ZTNB_HPP
#define ZTNB_HPP

#include <smithlab_utils.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>


class ZTNBD{
public:
  ZTNBD(const double m_, const double a_):
    mu(m_), alpha(a_) {};
  double get_mu() const {return mu;}
  double get_alpha() const {return alpha;}
  
  void set_mu(const double m_) {mu = m_;}
  void set_alpha(const double a_) {alpha = a_;}

  double operator()(int val) const;
  void estim_params(const std::vector<double> &vals_hist);
  void estim_params(const std::vector<double> &weighted_vals_hist,
                    const std::vector<double> &probs);
  double trunc_log_pdf(const size_t val);
  double log_pdf(const size_t val);

  double trunc_log_L(const std::vector<double> &vals_hist);

  double expected_zeros(const double distinct);
  double EM_estim_params(const double tol, const size_t max_iter,
			  std::vector<double> &vals_hist); //returns log_like
  double EM_estim_mu_fixed_alpha(const double tol, const size_t max_iter,
			       const std::vector<double> &vals_hist);
  double trunc_pval(const size_t val);
  // take 1- \sum_{X < val} P(X)
  double expected_inverse_sum(const double initial_distinct,
			      const double t);
  double expected_distinct(const double initial_distinct,
			   const double t);
  double expected_mincount(const size_t mincount,
			   const double initial_distinct,
			   const double t);
  //computes the expected # summands needed to reach a sum of "sum"

private:
  static const double max_allowed_alpha;
  static const double min_allowed_alpha;
  static const double tolerance;

  double score_fun_first_term(const std::vector<double> &pseudo_hist,
			      const double a_mid);
  double alpha_score_function(const std::vector<double> &pseudo_hist,
			      const double mean, const double a_mid,
			      const double vals_count);

  
  void set_helpers();
  
  double mu;
  double alpha;
  
  double n_helper;
  double p_helper;
  double n_log_p_minus_lngamma_n_helper;
  double log_q_helper;

};


#endif
