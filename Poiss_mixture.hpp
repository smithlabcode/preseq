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

#ifndef POISS_MIXTURE_HPP
#define POISS_MIXTURE_HPP


#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


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




class Poiss{
public:
  Poiss(const double l_): lambda(l_) {};
  double get_lambda() const {return(lambda);};
  double pdf(size_t val) const {
    return(gsl_ran_poisson_pdf(val, lambda));};

  double log_pdf(size_t val) const {
    return(-lambda + val*log(lambda) - gsl_sf_lnfact(val));};
  /*need to compute log_pdf directly, as gsl_ran_poisson
    rounds too quickly to 0 */


  void estim_param(const vector<size_t> &vals_hist, 
                   const vector<double> &probs);
  void estim_param(const vector<size_t> &vals_hist);
  string tostring() const {return toa(lambda);}
  void set_lambda(double lambda_hat) { lambda = lambda_hat;}
private:
  double lambda;

  static const double max_allowed_lambda;
  static const double min_allowed_lambda;
  static const double tolerance;

};

const double Poiss::max_allowed_lambda = 10000;
const double Poiss::min_allowed_lambda = 1e-20;
const double Poiss::tolerance = 1e-10; 

class ZTP : public Poiss
{
public:
  double trunc_pdf(size_t val) const {
    const double lamb = get_lambda();
    return(gsl_ran_poisson_pdf(val, lamb)/(1-exp(-lamb)));};

  double trunc_log_pdf(size_t val) const {
    const double lamb = get_lambda();
    return(-lamb + val*log(lamb) - gsl_sf_lnfact(val) -
	       log(1 - exp(-lamb)));};
  double trunc_log_L(const vector<size_t> &vals_hist);
  double EM_estim_param(const double tol, const size_t max_iter,
			vector<size_t> &vals_hist);
  double expected_zeros(const double pseudo_size);
  void trunc_estim_param_bisec(const double mean, const double tol);
  double expected_inverse_sum(const size_t sample_size,
			      const size_t sum);
  //computes the expected # summands needed to reach a sum of "sum"


};

class Poiss_mixture{
public:
  Poiss_mixture(const vector<Poiss> d_,
		const vector<double> m_): distros(d_), mixing(m_){;}

  vector<Poiss> get_distros() const {return distros;}
  vector<double> get_mixing() const {return mixing;}
  void set_distros(const vector<Poiss> dist) { distros = dist;}
  void set_mixing(const vector<double> mix) {mixing = mix;}
  double expectation_step(const vector<size_t> &vals_hist,
			  vector< vector<double> > &probs);
  void calculate_mixing(const vector<size_t> &vals_hist,
			const vector< vector<double> > &probs);
  void maximization_step(const vector<size_t> &vals_hist,
			 const vector< vector<double> > &probs);
  double log_L(const vector<size_t> &vals_hist);
  double EM_mix_resolve(const vector<size_t> &vals_hist,
			const double tol, const size_t max_iter);


private:
  vector<Poiss> distros;
  vector<double> mixing;
};

class ZTP_mixture{
public:
  ZTP_mixture(const vector<ZTP> d_,
	      const vector<double> mix_):
  distros(d_), mixing(mix_) {;}

  vector<ZTP> get_distros() const {return distros;}
  vector<double> get_mixing() const {return mixing;}
  void set_distros(const vector<ZTP> dist) { distros = dist;}
  void set_mixing(const vector<double> mix) {mixing = mix;}

  double trunc_expectation_step(const vector<size_t> &vals_hist,
				vector< vector<double> > &probs);
  void trunc_calculate_mixing(const vector<size_t> &vals_hist,
			const vector< vector<double> > &probs);
  void trunc_max_step(const vector<size_t> &vals_hist,
		      const vector< vector<double> > &probs,
		      const double tol);
  double trunc_log_L(const vector<size_t> &vals_hist);
  double EM_mix_resolve(const vector<size_t> &vals_hist,
			const double tol, const size_t max_iter);
  double expected_inverse_sum(const size_t sample_size,
			      const size_t sum);
  double Fisher_obs_info_2nd_theta(const vector< vector<double> > &probs,
				   const vector<size_t> &vals_hist,
				   const size_t indx);
  double Fisher_obs_info_mixed_theta(const vector< vector<double> > &probs,
				     const vector<size_t> &vals_hist,
				     const size_t indx1,
				     const size_t indx2);
  double Fisher_obs_info_2nd_lambda(const vector< vector<double> > &probs,
				    const vector<size_t> &vals_hist,
				    const size_t indx);
  double Fisher_obs_info_mixed_lambda(const vector< vector<double> > &probs,
				      const vector<size_t> &vals_hist,
				      const size_t indx1,
				      const size_t indx2);
  double Fisher_obs_info_mixed_same_indx(const vector< vector<double> > &probs,
					 const vector<size_t> &vals_hist,
					 const size_t indx);
  double Fisher_obs_info_mixed_last_indx(const vector< vector<double> > &probs,
				       const vector<size_t> &vals_hist,
				       const size_t theta_indx);
  double Fisher_obs_info_mixed_mixed_indx(const vector< vector<double> > &probs,
					  const vector<size_t> &vals_hist,
					  const size_t lambda_indx,
					  const size_t theta_indx);
  void compute_observed_Fisher_info(const vector<size_t> &vals_hist,
				    const vector< vector<double> > &probs);
private:
  vector<ZTP> distros;
  vector<double> mixing;
  vector< vector<double> > Fisher_info;
};

#endif
