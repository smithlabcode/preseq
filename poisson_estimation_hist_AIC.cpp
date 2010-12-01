#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

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


double 
log_sum_log_vec(const vector<double> &vals, size_t limit){
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
// abort if the sum is infinte //
#endif
    }
  }
  return(max_val + log(sum));
}

class Poiss{
public:
  Poiss(const double l_): lambda(l_) {};
  double get_lambda() const {return(lambda);};
  double get_pdf(size_t val) const {
    return(gsl_ran_poisson_pdf(val, lambda));};

  double get_log_pdf(size_t val) const {
    return(-lambda + val*log(lambda) - gsl_sf_lnfact(val));};
  /*need to compute log_pdf directly, as gsl_ran_poisson
    rounds too quickly to 0 */
  double get_nonzero_pdf(size_t val) const {
    return(gsl_ran_poisson_pdf(val, lambda)/(1-exp(-lambda)));};

  double get_nonzero_log_pdf(size_t val) const {
    return(-lambda + val*log(lambda) - gsl_sf_lnfact(val) -
	       log(1 - exp(-lambda)));};

  void estim_param(const vector<size_t> &vals_hist, 
                   const vector<double> &probs);
  void estim_param(const vector<size_t> &vals_hist);
  void estim_nonzero_param_bisec(const double mean);
  void estim_nonzero_param_newton(const double mean);
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


static double 
compute_mean(const vector<size_t> &vals_hist){
  double mean = 0.0;
  const double vals_size = static_cast<double>(
                       accumulate(vals_hist.begin(), 
                                  vals_hist.end(), 0));
  for(size_t i = 0; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}

static double
compute_weighted_mean(const vector<size_t> &values,
             const vector<double> &probs){
  return(std::inner_product(probs.begin(), probs.end(), 
			    values.begin(), 0.0)/
	 std::accumulate(probs.begin(), probs.end(), 0.0));
}

static double
compute_weighted_mean_hist(const vector<size_t> &vals_hist,
             const vector<double> &probs){
  double mean = 0.0;
  double probs_weight = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++){
    mean += i*probs[i]*vals_hist[i];
    probs_weight += vals_hist[i]*probs[i];
  }
  return(mean/probs_weight);
}

void 
Poiss::estim_param(const vector<size_t> &vals_hist){
  lambda = compute_mean(vals_hist);
  assert(finite(lambda));
}

void 
Poiss::estim_param(const vector<size_t> &vals_hist, 
                   const vector<double> &probs){
  lambda = compute_weighted_mean_hist(vals_hist,probs);

}
// lambda_k = weighted average(P(x_i \in  kth group)*x_i) //

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

static inline double
lambda_score_funct(const double mean, const double lambda){
  return(lambda + mean*exp(-lambda) - mean);
}

void
Poiss::estim_nonzero_param_bisec(const double mean){

  double lambda_low = mean-1;
  double lambda_high = mean;
  double lambda_mid = mean - 0.5;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while(movement(lambda_high, lambda_low) > tolerance &&
         diff > tolerance){
    
    lambda_mid = (lambda_low + lambda_high)/2;
    
    const double mid_val = lambda_score_funct(mean, lambda_mid);
    
    if (mid_val < 0) lambda_low = lambda_mid;
    else lambda_high = lambda_mid;
    
    diff = fabs((prev_val - mid_val)/max(mid_val, prev_val));
    
    prev_val = mid_val;
  }
  lambda = lambda_mid;
}

/* static double
lambda_newton_update(const double mean, const double lambda_hat){

  double deriv = 1 - mean*exp(-lambda_hat);
  return(lambda_hat - 
	 lambda_score_funct(mean,lambda_hat)/deriv);
}

void
Poiss::estim_nonzero_param_newton(const double mean){

  double lambda_hat = mean;
  double prev_lambda_hat = numeric_limits<double>::max();
  while(movement(lambda_hat, prev_lambda_hat) > tolerance){
    prev_lambda_hat = lambda_hat;
    lambda_hat = lambda_newton_update(mean, lambda_hat);
  }
  lambda = lambda_hat;
} */

static double
nonzero_expectation_step(const vector<size_t> &vals_hist,
                         const vector<Poiss> &distros,
                         const vector<double> &mixing,
                         vector< vector<double> > &probs){
  double score = 0.0;

  for(size_t i = 0; i < vals_hist.size(); i++){
    vector<double> log_denom_vec;
    
    for(size_t j = 0; j < distros.size(); j++){
      log_denom_vec.push_back(log(mixing[j]) +
                              distros[j].get_nonzero_log_pdf(i));
      if(!finite(log_denom_vec[j])){
        log_denom_vec[j] = -std::numeric_limits<double>::max()/
	  static_cast<double>(distros.size()); // arbitrary
      }
    }
    
    const double log_denom = log_sum_log_vec(log_denom_vec, 
                                             log_denom_vec.size());
    assert(finite(log_denom));
  
    for(size_t j = 0; j < distros.size(); j++)
      probs[j][i] = exp(log_denom_vec[j] - log_denom);

    score += i*log_denom;
  }
  return(score);
}

static void
calculate_mixing(const vector<size_t> &vals_hist,
                 const vector< vector<double> > &probs,
                 vector<double> &mixing){
  double vals_size = static_cast<double>(accumulate(vals_hist.begin(),
                                                    vals_hist.end(), 0));

  for(size_t i = 0; i < mixing.size(); i++){
    double inner_sum = 0.0;
    for(size_t j = 0; j < vals_hist.size(); j ++)
      inner_sum += vals_hist[j]*probs[i][j];
    mixing[i] = inner_sum/vals_size;
    assert(finite(mixing[i]));
  }
  double mixing_sum = accumulate(mixing.begin(), mixing.end(), 0.0);
  for(size_t i = 0; i < mixing.size(); i++){
    mixing[i] = mixing[i]/mixing_sum;
    assert(mixing[i] >= 0);
  }
}
    

static void
nonzero_maximization_step(const vector<size_t> &vals_hist,
                          const vector< vector<double> > &probs,
                          vector<double> &mixing,
                          vector<Poiss> &distros){

  for(size_t i = 0; i < distros.size(); i++){
    const double mean = compute_weighted_mean_hist(vals_hist, probs[i]);
    distros[i].estim_nonzero_param_bisec(mean);
  }

  calculate_mixing(vals_hist, probs, mixing);

}
// recalculate parameters, lambda and mixing_j = average(probs[j]) //

static double 
nonzero_logl_hist(const vector<size_t> &vals_hist,
                  const vector<Poiss> &distros,
                  const vector<double> &mixing){
  double logL = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++){
    double inner_sum = 0.0;
    vector<double> workspace_vec(distros.size(), 0.0);
    for(size_t j = 0; j < distros.size(); j++)
      workspace_vec[j] = log(mixing[j]) + 
                         distros[j].get_nonzero_log_pdf(i);
    
    inner_sum = log_sum_log_vec(workspace_vec, workspace_vec.size());
    logL += vals_hist[i]*inner_sum;
  }
  return(logL);
}

double 
nonzero_state_resolve(const vector<size_t> &vals_hist,
                      const double &tol,
                      const size_t max_iter,
                      vector<Poiss> &distros,
                      vector<double> &mixing){
  
  const size_t number_states = distros.size();
  double probs_starting_val = 1/static_cast<double>(number_states);
  vector< vector<double> > probs(number_states, 
                                 vector<double>(vals_hist.size(),
                                                probs_starting_val));

  double score = 0.0;
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_iter; ++i){
    cerr << "lambdas = ";
    for(size_t j = 0; j < distros.size(); j++)
      cerr << distros[j].get_lambda() << ", ";

    cerr << "\n mixings = ";
    for(size_t j = 0; j < mixing.size(); j++)
      cerr << mixing[j] << ", ";
    cerr << "\n";

    score = nonzero_expectation_step(vals_hist, distros, 
                                     mixing, probs);
    nonzero_maximization_step(vals_hist, probs, mixing, distros);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;
  }

  return(nonzero_logl_hist(vals_hist, distros, mixing));
}

static inline double
get_max(const double x1, const double x2){
  const double larger = (x1 > x2) ? x1 : x2;
  return(larger);
}

static inline double
calculate_neg_AIC(const size_t number_mixtures,
                  const double log_like){
  return(2*log_like - 2*(2*number_mixtures - 1));
}

static inline double
calculate_neg_BIC(const size_t number_mixtures,
                  const size_t sample_size,
                  const double log_like){
  return(2*log_like - (2*number_mixtures-1)*log(sample_size));
}

static inline size_t
select_number_mixtures(const vector<double> log_like){
  size_t current_max=0;
  double current_max_score = -std::numeric_limits<double>::max();
  for(size_t i = 0; i < log_like.size(); i++){
    const double holding_score = calculate_neg_AIC(i+1, log_like[i]);
    if(holding_score >= current_max_score){
      current_max_score = holding_score;
      current_max = i+1;
    }
  }
  return(current_max);
}
    			 

static void
get_counts(const vector<SimpleGenomicRegion> &read_locations,
	   vector<size_t> &values) {
  size_t count = 1;
  for (size_t i = 1; i < read_locations.size(); i++)
    if (read_locations[i] == read_locations[i-1])
      ++count;
    else {
      values.push_back(count);
      count = 1;
    }
} 
  


int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    size_t max_number_mixtures = 20;
    size_t max_iter = 100000;
    double tolerance = 1e-20;
 

    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("poisson_mix_estimation", "",
			   "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("max_number_mixtures", 'm', 
                      "number of mixtures to take",
                      false, max_number_mixtures);
    opt_parse.add_opt("tolerance", 't', "Numerical tolerance",
                     false, tolerance);


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

    vector<size_t> values;
    vector<SimpleGenomicRegion> read_locations;
    ReadBEDFile(input_file_name, read_locations);
    if (!check_sorted(read_locations))
      throw RMAPException("read_locations not sorted");

    get_counts(read_locations, values);
    size_t values_size = values.size();


    const size_t max_value = *std::max_element(values.begin(),values.end());
    vector<size_t> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i){
      ++vals_hist[static_cast<size_t>(values[i])];
    }

    vector<double> score_vec;
    vector< vector<double> > lambdas_vec;
    vector< vector<double> > mixings_vec;

    for(size_t k = 1; k <= max_number_mixtures; k++){
      size_t number_mixtures = k;
      vector<Poiss> distros;
      vector<double> mixing(number_mixtures, 
                          1/static_cast<double>(number_mixtures));
    
      sort(values.begin(), values.end());
      double step_size = floor(static_cast<double>(values.size())
                                 /static_cast<double>(number_mixtures));
      for(size_t j = 0; j < number_mixtures; j++){
        double lambda_hat = accumulate(values.begin()+j*step_size,
                              values.begin()+(j+1)*step_size, 
                              0.0)/step_size;

        if(j != 0){
          Poiss last_distro = distros.back();
          lambda_hat = get_max(lambda_hat, last_distro.get_lambda()) 
                         + 0.01;
        }
        else
          lambda_hat = lambda_hat+0.01;

 
        Poiss holding_distro(lambda_hat);
        distros.push_back(holding_distro);
      }
      double score = 0.0;  
      score  = nonzero_state_resolve(vals_hist, tolerance, max_iter, 
                                   distros, mixing);
      vector<double> lambdas;
      for(size_t j = 0; j < distros.size(); j++)
        lambdas.push_back(distros[j].get_lambda());

      lambdas_vec.push_back(lambdas);
      mixings_vec.push_back(mixing);
      score_vec.push_back(score);
    }

    size_t opt_num_mixs = select_number_mixtures(score_vec);
    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    *out << values.size() << endl;
    *out << opt_num_mixs << endl;
    *out << score_vec[opt_num_mixs-1] - (2*opt_num_mixs - 1) << endl;
    for(size_t j = 0; j < lambdas_vec[opt_num_mixs-1].size(); j++){
      *out << lambdas_vec[opt_num_mixs-1][j] << endl;
    } 

    for(size_t j = 0; j < mixings_vec[opt_num_mixs-1].size(); j++){
      *out << mixings_vec[opt_num_mixs-1][j] << endl;
    }
    
    

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

