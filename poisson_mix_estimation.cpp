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
  double get_pdf(unsigned int val) const {
    return(gsl_ran_poisson_pdf(val, lambda));};

  double get_log_pdf(unsigned int val) const {
    return(log(gsl_ran_poisson_pdf(val,lambda)));};

  double get_nonzero_pdf(unsigned int val) const {
    return(gsl_ran_poisson_pdf(val, lambda)/(1-exp(-lambda)));};

  double get_nonzero_log_pdf(unsigned int val) const {
    return(log(gsl_ran_poisson_pdf(val,lambda))-
	       log(1 - exp(-lambda)));};

  void estim_param(const vector<unsigned int> &values, 
                   const vector<double> &probs);
  void estim_param(const vector<unsigned int> &values);
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
compute_weighted_mean(const vector<unsigned int> &values,
             const vector<double> &probs){
  return(std::inner_product(probs.begin(), probs.end(), 
			    values.begin(), 0.0)/
	 std::accumulate(probs.begin(), probs.end(), 0.0));
}

void 
Poiss::estim_param(const vector<unsigned int> &values){
  lambda = accumulate(values.begin(), values.end(), 0.0)/
             static_cast<double>(values.size());
  assert(finite(lambda));
}

void 
Poiss::estim_param(const vector<unsigned int> &values, 
                   const vector<double> &probs){
  lambda = compute_weighted_mean(values,probs);

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

static double
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
}
    

static double
general_expectation_step(const vector<Poiss> &distros, 
                         const vector<double> &mixing,
                         const vector<unsigned int> &values, 
                         vector< vector<double> > &probs){
  double score = 0.0;

  for(size_t i = 0; i < values.size(); i++){ 
    vector<double> log_denom_vec;
  
    for(size_t j = 0; j < distros.size(); j++)
      log_denom_vec.push_back(log(mixing[j]) + 
                              distros[j].get_log_pdf(values[i]));

    double log_denom = log_sum_log_vec(log_denom_vec, log_denom_vec.size());

    assert(finite(log_denom));

    for(size_t j = 0; j < distros.size(); j++)
      probs[j][i] = exp(log_denom_vec[j] - log_denom);

    score += log_denom;

  }
  return(score);   
}
//calculate P(x_i in jth group) //


static void 
calculate_mixing(const vector<unsigned int> &values, 
                 const vector< vector<double> > &probs,
                 vector<double> &mixing){

  for(size_t i = 0; i < mixing.size(); i++){
    vector<double> log_probs;
    for(size_t j = 0; j < values.size(); j++){
      if(probs[i][j] != 0 ){
        log_probs.push_back(log(probs[i][j]));
      }
    }
    sort(log_probs.begin(), log_probs.end());
    mixing[i] = log_sum_log_vec(log_probs, log_probs.size()) - log(probs.size());
    assert(finite(mixing[i]));
  }
  double log_mixing_sum = log_sum_log_vec(mixing, mixing.size());
  for(size_t i = 0; i < mixing.size(); i++){
    mixing[i] = exp(mixing[i] - log_mixing_sum);
    assert(mixing[i] >= 0);
  }

}

static void
general_maximization_step(vector<unsigned int> &values,  
                          const vector< vector<double> > &probs, 
                          vector<double> &mixing, 
                          vector<Poiss> &distros){

  for(size_t i = 0; i < distros.size(); i++)
    distros[i].estim_param(values, probs[i]);

  calculate_mixing(values, probs, mixing);

}
// recalculate parameters, lambda and mixing_j = average(probs[j]) //

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

double
general_state_resolve(const double &tol, 
                      const size_t &max_itr,
                      vector<unsigned int> &values, 
                      vector<Poiss> &distros, 
                      vector<double> &mixing){


  const size_t number_states = distros.size();
  double probs_starting_val = 1/static_cast<double>(number_states);
  vector< vector<double> > probs(number_states, vector<double>(values.size(), 
                                 probs_starting_val));

  double score = 0.0;
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_itr; ++i){
    cerr << "lambdas = ";
    for(size_t j = 0; j < distros.size(); j++)
      cerr << distros[j].get_lambda() << ", ";

    cerr << "\n mixings = ";
    for(size_t j = 0; j < mixing.size(); j++)
      cerr << mixing[j] << ", ";
    cerr << "\n";
    score = general_expectation_step(distros, mixing, values, probs);
    general_maximization_step(values, probs, mixing, distros);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;
    cerr << "lambdas = ";
    for(size_t j = 0; j < distros.size(); j++)
      cerr << distros[j].get_lambda() << ", ";

    cerr << "\n mixings = ";
    for(size_t j = 0; j < mixing.size(); j++)
      cerr << mixing[j] << ", ";
    cerr << "\n";
  }

  return(general_log_l(values, distros, mixing));
}

static double
general_nonzero_expectation_step(const vector<Poiss> &distros, 
                                 const vector<double> &mixing,
                                 const vector<unsigned int> &values, 
                                 vector< vector<double> > &probs){
  double score = 0.0;

  for(size_t i = 0; i < values.size(); i++){ 
    vector<double> log_denom_vec;
  
    for(size_t j = 0; j < distros.size(); j++)
      log_denom_vec.push_back(log(mixing[j]) + 
                              distros[j].get_nonzero_log_pdf(values[i]));

    double log_denom = log_sum_log_vec(log_denom_vec, log_denom_vec.size());

    assert(finite(log_denom));

    for(size_t j = 0; j < distros.size(); j++)
      probs[j][i] = exp(log_denom_vec[j] - log_denom);

    score += log_denom;

  }
  return(score);   
}

static void
general_nonzero_maximization_step(vector<unsigned int> &values,  
                                  const vector< vector<double> > &probs, 
                                  vector<double> &mixing, 
                                  vector<Poiss> &distros){

  for(size_t i = 0; i < distros.size(); i++){
    double mean = compute_weighted_mean(values, probs[i]);
    distros[i].estim_nonzero_param_bisec(mean);
  }

  calculate_mixing(values, probs, mixing);

}
// recalculate parameters, lambda and mixing_j = average(probs[j]) //

static double 
general_nonzero_log_l(const vector<unsigned int> &values, 
                      const vector<Poiss> &distros,
                      const vector<double> &mixing){

  vector<double> score_vec;

  for(size_t i = 0; i < values.size(); i++){
    double inner_sum = 0.0;
    vector<double> workspace_sum(distros.size(), 0.0);

    for(size_t j = 0; j < mixing.size(); j++)
      workspace_sum[j] = log(mixing[j]) + 
                           distros[j].get_nonzero_log_pdf(values[i]);

    inner_sum = log_sum_log_vec(workspace_sum, workspace_sum.size());
    score_vec.push_back(inner_sum);
  }
  sort(score_vec.begin(), score_vec.end());
  return(accumulate(score_vec.begin(), score_vec.end(), 0.0));
}


static double
complete_nonzero_log_l(const vector<unsigned int> &values, 
                       const vector<Poiss> &distros,
                       const vector< vector<double> > &probs){
  
  vector<double> score_vec;

  for(size_t i = 0; i < values.size(); i ++){
    vector<double> log_score;

    for(size_t j = 0; j < distros.size(); j++){
      if(probs[j][i] > 0){
        log_score.push_back(log(probs[j][i]) + 
                            distros[i].get_nonzero_log_pdf(values[i]));
      }
    }
   
    score_vec.push_back(log_sum_log_vec(log_score, log_score.size()));
  }
  
  return(accumulate(score_vec.begin(), score_vec.end(), 0.0));
}

double
general_nonzero_state_resolve(const double &tol, 
                              const size_t &max_itr,
                              vector<unsigned int> &values, 
                              vector<Poiss> &distros, 
                              vector<double> &mixing){


  const size_t number_states = distros.size();
  double probs_starting_val = 1/static_cast<double>(number_states);
  vector< vector<double> > probs(number_states, vector<double>(values.size(), 
                                 probs_starting_val));

  double score = 0.0;
  double error = std::numeric_limits<double>::max();
  double prev_score = std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_itr; ++i){
    cerr << "lambdas = ";
    for(size_t j = 0; j < distros.size(); j++)
      cerr << distros[j].get_lambda() << ", ";

    cerr << "\n mixings = ";
    for(size_t j = 0; j < mixing.size(); j++)
      cerr << mixing[j] << ", ";
    cerr << "\n";
    score = general_nonzero_expectation_step(distros, mixing, values, probs);
    general_nonzero_maximization_step(values, probs, mixing, distros);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;
  }

  return(general_nonzero_log_l(values, distros, mixing));
}


int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    size_t number_mixtures = 1;
    bool zeros = false;
    size_t max_itr = 100000;
    double tolerance = 1e-20;
    string lambdas_string;

    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("poisson_mix_estimation", "",
			   "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("number_mixtures", 'm', "number of mixtures to take",
                      false, number_mixtures);
    opt_parse.add_opt("tolerance", 't', "Numerical tolerance",
                     false, tolerance);
    opt_parse.add_opt("initial_lambdas", 'l', "initial guesses on lambdas",
                      false, lambdas_string);
    opt_parse.add_opt("zeros", 'z', "Include if zeros are included",
                      false, zeros);

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

    vector<unsigned int> values;
    std::ifstream in(input_file_name.c_str());
    if (!in) 
      throw BEDFileException("cannot open input file " + input_file_name);
    static const size_t buffer_size = 10000; // Magic!
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    double real_score = atof(buffer);
    while (!in.eof()){
      int holding_val = 0;
      char buffer[buffer_size];
      in.getline(buffer, buffer_size);
      holding_val = atoi(buffer);
      if(!(holding_val >= 0)){
        cerr << "Invalid Input.\n";
      }
      assert(holding_val >= 0);
      values.push_back(holding_val);
    }
    in.close();
    
    vector<Poiss> distros;
    vector<double> mixing(number_mixtures, 
                          1/static_cast<double>(number_mixtures));
    
    if(lambdas_string.size() == 0){
      sort(values.begin(), values.end());
      double step_size = floor(static_cast<double>(values.size())
                                 /static_cast<double>(number_mixtures));
      for(size_t i = 0; i < number_mixtures; i++){
        double lambda_hat = accumulate(values.begin()+i*step_size,
                              values.begin()+(i+1)*step_size, 0.0)/step_size;
        lambda_hat = lambda_hat + 0.01;
        Poiss holding_distro(lambda_hat);
        distros.push_back(holding_distro);
      }
               
                
    }
    else{
      vector<double> lambdas;
      vector<string> lambdas_parts = rmap::split(lambdas_string, ",");    
      for(size_t i = 0; i < number_mixtures; i++){
        lambdas.push_back(atof(lambdas_parts[i].c_str()));
        Poiss holding_distro(lambdas[i]);
        distros.push_back(holding_distro);
      }
    }

    double score = 0.0;

    if(zeros == true)
      score = general_state_resolve(tolerance, max_itr, values, distros, mixing);
    
    else
      score = general_nonzero_state_resolve(tolerance, max_itr, values,
                                            distros, mixing);
    

    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    *out << "lambdas" << endl;
    for(size_t j = 0; j < distros.size(); j++){
      *out << distros[j].get_lambda() << endl;
    } 
    *out << "mixings" << endl;
    for(size_t j = 0; j < mixing.size(); j++){
      *out << mixing[j] << endl;
    }
    *out << "estimated log likelihood" << "\t" << score << endl;
    *out << "real log likelihood" << "\t" << real_score << endl;
    
    

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


