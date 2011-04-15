/*    NBD_mixture:
 *
 *    Copyright (C) 2011 University of Southern California and
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


#ifndef NBD_MIXTURE_HPP
#define NBD_MIXTURE_HPP

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

using std::max;
using std::setw;
using std::fabs;
using std::ceil;
using std::greater;
using std::numeric_limits;
using std::cin;


class NBD {
public:
  NBD(const double m_, const double a_) : mu(m_), alpha(a_) {set_helpers();}

  double get_mu() const {return mu;}
  double get_alpha() const {return alpha;}
  
  void set_mu(const double m_) {mu = m_;}
  void set_alpha(const double a_) {alpha = a_;}
  
  double operator()(int val) const;
  void estim_params(const vector<size_t> &val_hist);
  void estim_params(const vector<size_t> &vals_hist,
		    const vector<double> &probs);
  double log_pdf(const size_t val);
  double log_L(const vector<size_t> &vals_hist);


 

private:
  
  static const double max_allowed_alpha;
  static const double min_allowed_alpha;
  static const double tolerance;

  double score_fun_first_term(const vector<size_t> &vals_hist,
			      const double a_mid);
  double alpha_score_function(const vector<size_t> &vals_hist,
			      const double mean, const double a_mid,
			      const double vals_count);
  double score_fun_first_term(const vector<double> &pseudo_hist,
			      const double a_mid);
  double alpha_score_function(const vector<double> &pseudo_hist,
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



class ZTNBD{
public:
  ZTNBD(const double m_, const double a_):
    mu(m_), alpha(a_) {};
  double get_mu() const {return mu;}
  double get_alpha() const {return alpha;}
  
  void set_mu(const double m_) {mu = m_;}
  void set_alpha(const double a_) {alpha = a_;}

  double operator()(int val) const;
  void estim_params(const vector<size_t> &val_hist);
  void estim_params(const vector<size_t> &vals_hist,
		    const vector<double> &probs);
  
  double trunc_log_pdf(size_t val);
  double trunc_log_L(const vector<size_t> &vals_hist);
  double expected_zeros(const double pseudo_size);
  double EM_estim_params(const double tol, const size_t max_iter,
			  vector<size_t> &vals_hist); //returns log_like
  double trunc_pval(const size_t val, const double tol);
  // tol=tolerance, take \sum_{X>=val} P(X), stop when diff < tol
  double expected_inverse_sum(const size_t sample_size,
			      const size_t sum);
  //computes the expected # summands needed to reach a sum of "sum"

private:
  static const double max_allowed_alpha;
  static const double min_allowed_alpha;
  static const double tolerance;

  double score_fun_first_term(const vector<size_t> &vals_hist,
			      const double a_mid);
  double alpha_score_function(const vector<size_t> &vals_hist,
			      const double mean, const double a_mid,
			      const double vals_count);
  double score_fun_first_term(const vector<double> &pseudo_hist,
			      const double a_mid);
  double alpha_score_function(const vector<double> &pseudo_hist,
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




class NBD_mixture{
public:
  NBD_mixture(const vector<NBD> d_, 
	      const vector<double> mix_): 
    distros(d_), mixing(mix_){;}

  vector<NBD> get_distros() const {return distros;}
  vector<double> get_mixing() const {return mixing;}
  void set_distros(const vector<NBD> dist) {distros = dist;}
  void set_mixing(const vector<double> mix) {mixing = mix;}
  void calculate_mixing(const vector<size_t> &vals_hist,
			const vector< vector<double> > &probs);
  void expectation_step(const vector<size_t> &vals_hist,
			vector< vector<double> > &probs);
  void maximization_step(const vector<size_t> &vals_hist,
			 const vector< vector<double> > &probs);
  double log_L(const vector<size_t> &vals_hist);
  double EM_resolve_mix(const vector<size_t> &vals_hist,
			const double &tol, const size_t max_iter);
private:
  vector<NBD> distros;
  vector<double> mixing;
};

class ZTNBD_mixture{
public:
  ZTNBD_mixture(const vector<ZTNBD> d_,
		const vector<double> mix_,
		const vector< vector<double> > Fish_):
    distros(d_), mixing(mix_), Fisher_info(Fish_){;}

  vector<ZTNBD> get_distros() const {return distros;}
  vector<double> get_mixing() const {return mixing;}
  vector< vector<double> > get_Fish_info() const {return Fisher_info;}
  void set_distros(const vector<ZTNBD> dist) {distros = dist;}
  void set_mixing(const vector<double> mix) {mixing = mix;}
  void set_Fish_info(const vector< vector<double> > Fish) {Fisher_info=Fish;}

  void trunc_calculate_mixing(const vector<size_t> &vals_hist,
			      const vector< vector<double> > &probs);
  void trunc_expectation_step(const vector<size_t> &vals_hist,
			      vector< vector<double> > &probs);
  void trunc_maximization_step(const vector<size_t> &vals_hist,
			       const vector< vector<double> > &probs);

  double trunc_log_L(const vector<size_t> &vals_hist);
  void compute_mixing_w_zeros(const size_t sample_size,
			      vector<double> &mixing);
  double EM_resolve_mix_add_zeros(const double &tol, const size_t max_iter,
				  vector<size_t> &vals_hist);
  double expected_inverse_sum(const size_t sample_size,
			      const size_t sum);
  void compute_Fisher_info(const vector< vector<double> > &probs,
			   const vector<size_t> &vals_hist,
			   const vector<double> &mixing,
			   const double expected_zeros);
private:
  vector<ZTNBD> distros;
  vector<double> mixing;
  vector< vector<double> > Fisher_info;
};

double 
log_sum_log_vec(const vector<double> &vals, size_t limit);

#endif
