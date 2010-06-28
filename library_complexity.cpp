/*    library_complexity: a program for estimating the complexity of a
 *    DNA library prepared for sequencing based on sequencing a small
 *    number of reads
 *
 *    Copyright (C) 2010 University of Southern California and
 *                       Andrew D. Smith
 *                       Timothy Dailey
 *
 *    Authors: Andrew D. Smith and Timothy Dailey
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

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include <fstream>
#include <iomanip>
#include <numeric>

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

double
log_sum_log_vec(const vector<double> &vals, size_t limit) {
  const size_t max_idx = 
    max_element(vals.begin(), vals.begin() + limit) - 
    vals.begin();
  const double max_val = vals[max_idx];
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(finite(sum));
#endif
    }
  }
  return max_val + log(sum);
}

class NBD {
public:
  NBD(const double m_, const double a_) : mu(m_), alpha(a_) {set_helpers();}

  double get_mu() const {return mu;}
  double get_alpha() const {return alpha;}
  
  void set_mu(const double m_) {mu = m_;}
  void set_alpha(const double a_) {alpha = a_;}
  
  double operator()(double val) const;
  void estim_params(const vector<double> &vals);
  void estim_params(const vector<double> &vals,
		    const vector<double> &probs);
  string tostring() const {return toa(mu) + "," + toa(alpha);}
private:
  
  static const double max_allowed_alpha;
  static const double min_allowed_alpha;
  static const double tolerance;
  
  void set_helpers();
  
  double mu;
  double alpha;
  
  double n_helper;
  double p_helper;
  double n_log_p_minus_lngamma_n_helper;
  double log_q_helper;
};

const double NBD::max_allowed_alpha = 100;
const double NBD::min_allowed_alpha = 1e-20;
const double NBD::tolerance = 1e-10; 

void
NBD::set_helpers() {
  n_helper = 1/alpha;
  p_helper = n_helper/(n_helper + mu);
  n_log_p_minus_lngamma_n_helper = 
    n_helper*log(p_helper) - gsl_sf_lngamma(n_helper);
  log_q_helper = log(1 - p_helper);
  //TODO: check these!!!!!!!!!!!!!!!!
}


// THIS IS THE LOG OF THE LIKELIHOOD!!!!!
double 
NBD::operator()(const double val) const {
  const double P = (gsl_sf_lngamma(val + n_helper) - 
		    gsl_sf_lnfact(static_cast<size_t>(val))) +
    n_log_p_minus_lngamma_n_helper + val*log_q_helper;
  if (!finite(P)) return -40;
  return P;
}

static inline double
score_fun_first_term(const vector<double> &vals_hist, 
		     const double mu, const double alpha) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
	inner_sum += j/(1 + alpha*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum;
}

static inline double
alpha_score_function(const vector<double> &vals_hist, const double mu, 
		     const double alpha, const double vals_count) {
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_fun_first_term(vals_hist, mu, alpha)/vals_count + 
	  (log(one_plus_alpha_mu)/alpha - mu)/alpha);
}

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b));
}

void
NBD::estim_params(const vector<double> &vals) {
  // This is the mu
  mu = accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
  
  // Now for the alpha
  const double max_value = *max_element(vals.begin(), vals.end());
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < vals.size(); ++i)
    ++vals_hist[static_cast<size_t>(vals[i])];
  
  const double vals_count = vals.size();
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;
  set_helpers();
}

void
NBD::estim_params(const vector<double> &vals, const vector<double> &probs) {
  const size_t lim = vals.size();

  vector<double> workspace_vals(lim), workspace_probs(lim);
  for (size_t i = 0; i < lim; ++i) {
    workspace_probs[i] = log(probs[i]);
    workspace_vals[i] = log(vals[i]) + log(probs[i]);
  }
  const double vals_count = exp(log_sum_log_vec(workspace_probs, lim));
  mu = exp(log_sum_log_vec(workspace_vals, lim))/vals_count;
  
  // now the alpha
  const double max_value = *max_element(vals.begin(), vals.begin() + lim);
  vector<double> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);
  for (size_t i = 0; i < lim; ++i)
    vals_hist[static_cast<size_t>(vals[i])] += probs[i];
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    
    a_mid = (a_low + a_high)/2;
    
    const double mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    
    diff = fabs((prev_val - mid_val)/max(mid_val, prev_val));
    
    prev_val = mid_val;
  }
  alpha = a_mid;
  set_helpers();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static double
expectation_step(const vector<double> &values, const double &mixing,
		 const NBD &distro1, const NBD &distro2,
		 vector<double> &probs1, vector<double> &probs2) {

  const double log_mixing1 = log(mixing);
  const double log_mixing2 = log(1 - mixing);
  
  assert(finite(log_mixing1));
  assert(finite(log_mixing2));
  
  double score = 0;
  vector<double>::iterator p1(probs1.begin()), p2(probs2.begin());
  for (vector<double>::const_iterator i(values.begin()); i != values.end(); ++i) {
    const double term1 = log_mixing1 + distro1(*i);
    const double term2 = log_mixing2 + distro2(*i);
    const double denom = ((term1 > term2) ?
			  term1 + log(1.0 + exp(term2 - term1)) :
			  term2 + log(1.0 + exp(term1 - term2)));

    assert(finite(term1));
    assert(finite(term2));
    assert(finite(denom));
    
    *p1++ = exp(term1 - denom);
    *p2++ = exp(term2 - denom);
    
    score += denom;
  }
  return score;
}

static void
maximization_step(const vector<double> &values, 
		  const vector<double> &probs1, const vector<double> &probs2,
		  double &mixing, NBD &distro1, NBD &distro2) {
  distro1.estim_params(values, probs1);
  distro2.estim_params(values, probs2);
  
  vector<double> log_probs1(probs1);
  vector<double> log_probs2(probs2);
  for (size_t i = 0; i < log_probs1.size(); ++i) {
    log_probs1[i] = log(log_probs1[i]);
    log_probs2[i] = log(log_probs2[i]);
  }    
  mixing = rmap::log_sum_log_vec(log_probs1, log_probs1.size());
  const double bg_mixing = rmap::log_sum_log_vec(log_probs2, log_probs2.size());
  const double mix_sum = ((mixing > bg_mixing) ?
			  mixing + log(1 + exp(bg_mixing - mixing)) :
			  bg_mixing + log(1 + exp(mixing - bg_mixing)));
  mixing = exp(mixing - mix_sum);
}


static string
info_line(const size_t itr, const double score, const double prev_score, 
	  const double mixing, NBD &distro1, NBD &distro2) {
  std::ostringstream ss;
  ss << itr << "\t"
     << setw(10) << std::setprecision(4) 
     << (prev_score - score)/prev_score << "\t"
     << setw(14) << distro1.tostring() << "\t" 
     << setw(10) << mixing << "\t"
     << setw(14) << distro2.tostring() << "\t" 
     << setw(10) << 1 - mixing;
  return ss.str();
}


void
TwoStateResolveMixture(const bool VERBOSE,
		       const size_t max_itr, const double tolerance, 
		       const vector<double> &values,
		       NBD &distro1, NBD &distro2, double &mixing) {
  
  vector<double> probs1(values.size(), 0.0), probs2(values.size(), 0.0);
  if (VERBOSE)
    cerr << setw(10) << "DELTA"
	 << setw(14) << "(PARAMS,MIX)" << endl;
  
  // expectation maximization
  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_itr; ++i) {
    const double score = 
      expectation_step(values, mixing, distro1, distro2, probs1, probs2);
    maximization_step(values, probs1, probs2, mixing, distro1, distro2);
    if (VERBOSE)
      cerr << info_line(i, score, prev_score, mixing, distro1, distro2) << endl;
    if ((prev_score - score)/prev_score < tolerance)
      break;
    prev_score = score;
  }
}


static void
get_counts(const vector<SimpleGenomicRegion> &read_locations,
	   vector<double> &values) {
  size_t count = 1;
  for (size_t i = 1; i < read_locations.size(); ++i)
    if (read_locations[i] != read_locations[i-1])
      ++count;
    else {
      values.push_back(count);
      count = 0;
    }
  values.push_back(count);
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    size_t max_itr = 250;
    double tolerance = 1e-20;
    double mixing = 0.5;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("library_complexity", 
			   "a program for estimating the complexity of a "
			   "DNA library prepared for sequencing based on "
			   "sequencing a small number of reads",
			   "<mapped-read-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("tolerance", 't', "numerical tolerance", 
		      false , tolerance);
    opt_parse.add_opt("itr", 'i', "max number of iterations", 
		      false , max_itr);
    opt_parse.add_opt("mixing", 'm', "initial guess for mixing parameter", 
		      false , mixing);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file_name = leftover_args.front();
    /**********************************************************************/

    vector<SimpleGenomicRegion> read_locations;
    ReadBEDFile(input_file_name, read_locations);
    if (!check_sorted(read_locations))
      throw RMAPException("read_locations not sorted");

    vector<double> values;
    get_counts(read_locations, values);
    
    NBD distro(0.0, 0.0);
    distro.estim_params(values);
    
    if (VERBOSE) 
      cerr << "FULL DATA SET:\t" << distro.tostring() << endl;
    
    NBD distro1(2*distro.get_mu(), distro.get_alpha());
    NBD distro2(distro.get_mu()/2, distro.get_alpha());
    
    TwoStateResolveMixture(VERBOSE, max_itr, tolerance, values,
			   distro1, distro2, mixing);
    
    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    *out << "DISTRO1:\t" << distro1.tostring() << endl
	 << "DISTRO2:\t" << distro2.tostring() << endl;
    if (out != &std::cout) delete out;
  }
  catch (RMAPException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
