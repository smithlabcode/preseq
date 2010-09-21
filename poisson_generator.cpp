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
using std::cin;


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
  double get_pdf(unsigned int val) const {
    return(gsl_ran_poisson_pdf(val, lambda));};

  double get_log_pdf(unsigned int val) const {
    return(val*log(lambda) - lambda - gsl_sf_lnfact(val));};

  double get_trun_pdf(unsigned int val) const {
    return(gsl_ran_poisson_pdf(val, lambda)/(1-exp(-lambda)));};

  double get_trun_log_pdf(unsigned int val) const {
    return(val*log(lambda) - lambda - gsl_sf_lnfact(val) 
           - log(1-exp(-lambda)));};

  void estim_param(const vector<unsigned int> &values, const vector<double> &probs);
  void estim_param(const vector<unsigned int> &values);
  void estim_trun_param(const vector<unsigned int> &values, const vector<double> &probs);
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

void 
Poiss::estim_param(const vector<unsigned int> &values){
  lambda = accumulate(values.begin(), values.end(), 0.0)/
             static_cast<double>(values.size());
  assert(finite(lambda));
}

void 
Poiss::estim_param(const vector<unsigned int> &values, 
                   const vector<double> &probs){
  lambda =  std::inner_product(probs.begin(), probs.end(), 
			    values.begin(), 0.0)/
    std::accumulate(probs.begin(), probs.end(), 0.0);

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
Poiss::estim_trun_param(const vector<unsigned int> &values,
                             const vector<double> &probs){
  
  double mean = std::inner_product(probs.begin(), probs.end(), 
			           values.begin(), 0.0)/
                std::accumulate(probs.begin(), probs.end(), 0.0);

  double lambda_low = min_allowed_lambda;
  double lambda_high = max_allowed_lambda;
  double lambda_mid = max_allowed_lambda;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(lambda_high, lambda_low) > tolerance) {
    
    lambda_mid = (lambda_low + lambda_high)/2;
    
    const double mid_val = lambda_score_funct(mean, lambda_mid);
    
    if (mid_val < 0) lambda_low = lambda_mid;
    else lambda_high = lambda_mid;
    
    diff = fabs((prev_val - mid_val)/max(mid_val, prev_val));
    
    prev_val = mid_val;
  }
  lambda = lambda_mid;
}

void
generate_mix_poisson(const vector<double> &mixing,
                     const vector<Poiss> &distros,
                     vector<unsigned int> &sample){

  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());

  for(size_t i = 0; i < sample.size(); i++){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < mixing.size(); j++){
      test_val += mixing[j];
      if(u <= test_val){
        sample[i] = gsl_ran_poisson(rng, distros[j].get_lambda());
        break;
      }
    }
  }
     
} 

void
generate_mix_trun_poisson(const vector<double> &mixing,
                     const vector<Poiss> &distros,
                     vector<unsigned int> &sample){

  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());

  for(size_t i = 0; i < sample.size(); i++){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < mixing.size(); j++){
      test_val += mixing[j];
      if(u <= test_val){
        sample[i] = gsl_ran_poisson(rng, distros[j].get_lambda());
        while(sample[i] == 0){
          sample[i] = gsl_ran_poisson(rng, distros[j].get_lambda());
        }
        break;
      }
    }
  }

}     

static double 
general_log_l(const vector<unsigned int> &values, 
              const vector<Poiss> &distros,
              const vector<double> &mixing){

  vector<double> score_vec;

  for(size_t i = 0; i < values.size(); i++){
    double inner_sum = 0.0;
    vector<double> workspace_sum(distros.size(), 0.0);

    for(size_t j = 0; j < mixing.size(); j++)
      workspace_sum[j] = log(mixing[j]) + distros[j].get_log_pdf(values[i]);

    inner_sum = log_sum_log_vec(workspace_sum, workspace_sum.size());
    score_vec.push_back(inner_sum);
  }
  sort(score_vec.begin(), score_vec.end());
  return(accumulate(score_vec.begin(), score_vec.end(), 0.0));
}


static double
complete_log_l(const vector<unsigned int> &values, 
              const vector<Poiss> &distros,
              const vector< vector<double> > &probs){
  
  vector<double> score_vec;

  for(size_t i = 0; i < values.size(); i ++){
    vector<double> log_score;

    for(size_t j = 0; j < distros.size(); j++){
      if(probs[j][i] > 0){
        log_score.push_back(log(probs[j][i]) + 
                            distros[i].get_log_pdf(values[i]));
      }
    }
   
    score_vec.push_back(log_sum_log_vec(log_score, log_score.size()));
  }
  
  return(accumulate(score_vec.begin(), score_vec.end(), 0.0));
}

static double 
general_trun_log_l(const vector<unsigned int> &values, 
              const vector<Poiss> &distros,
              const vector<double> &mixing){

  vector<double> score_vec;

  for(size_t i = 0; i < values.size(); i++){
    double inner_sum = 0.0;
    vector<double> workspace_sum(distros.size(), 0.0);

    for(size_t j = 0; j < mixing.size(); j++)
      workspace_sum[j] = log(mixing[j]) + distros[j].get_trun_log_pdf(values[i]);

    inner_sum = log_sum_log_vec(workspace_sum, workspace_sum.size());
    score_vec.push_back(inner_sum);
  }
  sort(score_vec.begin(), score_vec.end());
  return(accumulate(score_vec.begin(), score_vec.end(), 0.0));
}


static double
complete_trun_log_l(const vector<unsigned int> &values, 
              const vector<Poiss> &distros,
              const vector< vector<double> > &probs){
  
  vector<double> score_vec;

  for(size_t i = 0; i < values.size(); i ++){
    vector<double> log_score;

    for(size_t j = 0; j < distros.size(); j++){
      if(probs[j][i] > 0){
        log_score.push_back(log(probs[j][i]) + 
                            distros[i].get_trun_log_pdf(values[i]));
      }
    }
   
    score_vec.push_back(log_sum_log_vec(log_score, log_score.size()));
  }
  
  return(accumulate(score_vec.begin(), score_vec.end(), 0.0));
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    double sample_size = 500000;
    size_t number_mixtures = 1;
    bool zeros = false;
    string lambdas_string;
    string mixing_string;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("poisson_generator", "",
			   "");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("sample_size", 'n', "sample size to take", 
		      false , sample_size);
    opt_parse.add_opt("zeros", 'z', "zeros are included",
                      false, zeros);
    opt_parse.add_opt("lambdas", 'l', "lambdas, separated by commas",
                      true, lambdas_string);
   opt_parse.add_opt("mixing", 'x', "mixing parameters, separated by commas",
                      true, mixing_string);
    opt_parse.add_opt("number_mixtures", 'm', "number of mixtures",
                      false, number_mixtures);
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


    /**********************************************************************/

    vector<Poiss> distros;
    vector<double> mixing;
    vector<double> lambdas;

    vector<string> lambdas_parts = rmap::split(lambdas_string, ",");    
                                                      
    vector<string> mixing_parts = rmap::split(mixing_string, ",");



    for(size_t i = 0; i < number_mixtures; i++){
      lambdas.push_back(atof(lambdas_parts[i].c_str()));
      mixing.push_back(atof(mixing_parts[i].c_str()));
    }

    if(accumulate(mixing.begin(), mixing.end(), 0.0) != 1){
      cerr << "Mixing parameters must sum to 1! \n";
    }
    assert(accumulate(mixing.begin(), mixing.end(), 0.0) == 1);

    for(size_t i = 0; i < lambdas.size(); i++){
      Poiss holding_distro(lambdas[i]);
      distros.push_back(holding_distro);
    }

    vector<unsigned int> the_sample(sample_size, 0);
    double score = 0.0;


    if(zeros == true){
      generate_mix_poisson(mixing, distros, the_sample);

      score = general_log_l(the_sample, distros, mixing);
    }
    else{
      generate_mix_trun_poisson(mixing, distros, the_sample);

      score = general_trun_log_l(the_sample, distros, mixing);
    }

    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    *out << score << endl;
    for(size_t j = 0; j < the_sample.size(); j++){
      *out << the_sample[j] << endl;
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
      


