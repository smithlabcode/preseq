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




int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    string combine_string;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("extrapolate_library_complexity", "",
			   "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("combine_which", 'c', "which components to combine",
                      true, combine_string);

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

    vector<string> combine_parts = rmap::split(combine_string, ",");
    vector<size_t> combine_vec;
    for(size_t i = 0; i < combine_parts.size(); i++)
      combine_vec.push_back(atoi(combine_parts[i].c_str()));



    vector<double> lambdas_vec;
    vector<double> mixing_vec;

    std::ifstream in(input_file_name.c_str());
    if (!in) 
      throw BEDFileException("cannot open input file " 
                               + input_file_name);

    static const size_t buffer_size = 10000; // Magic!
    char buffer1[buffer_size];
    in.getline(buffer1, buffer_size);
    size_t values_size = atoi(buffer1);
    char buffer2[buffer_size];
    in.getline(buffer2, buffer_size);
    size_t number_mixtures = atoi(buffer2);
    char buffer3[buffer_size];
    in.getline(buffer3, buffer_size);
    double log_like = atof(buffer3);
    size_t indx = 0;
    while (!in.eof() && indx < number_mixtures){
      double holding_val = 0;
      char buffer[buffer_size];
      in.getline(buffer, buffer_size);
      holding_val = atof(buffer);
      if(!(holding_val > 0)){
        cerr << "Invalid Input.\n";
      }
      assert(holding_val > 0);
      lambdas_vec.push_back(holding_val);
      indx++;
    }
    indx = 0;
    while (!in.eof() && indx < number_mixtures){
      double holding_val = 0;
      char buffer[buffer_size];
      in.getline(buffer, buffer_size);
      holding_val = atof(buffer);
      if(!(holding_val > 0)){
        cerr << "Invalid Input.\n";
      }
      assert(holding_val > 0);
      mixing_vec.push_back(holding_val);
      indx++;
    }
    in.close();

    double mean = 0.0;
    double mixing = 0.0;
    for(size_t i = 0; i < combine_vec.size(); i++){
      mixing += mixing_vec[combine_vec[i]-1];
      mean += mixing_vec[combine_vec[i]-1]*lambdas_vec[combine_vec[i]-1]/
	(1 - exp(- lambdas_vec[combine_vec[i]-1]));
    }
    mean = mean/mixing;
    Poiss distro(2);
    distro.estim_nonzero_param_bisec(mean);

    for(size_t i = combine_vec.size()-1; i > 0 ; --i){
      mixing_vec.erase(mixing_vec.begin()+combine_vec[i]-1);
      lambdas_vec.erase(lambdas_vec.begin()+combine_vec[i] - 1);
    }
    mixing_vec.erase(mixing_vec.begin()+combine_vec[0]-1);
    lambdas_vec.erase(lambdas_vec.begin()+combine_vec[0]-1);
    
    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    *out << values_size << endl;
    *out << number_mixtures - combine_vec.size() + 1 << endl;
    *out << log_like << endl;
    *out << distro.get_lambda() << endl;
    for(size_t i = 0; i < lambdas_vec.size(); i++)
      *out << lambdas_vec[i] << endl;
    *out << mixing << endl;
    for(size_t i = 0; i < mixing_vec.size(); i++)
      *out << mixing_vec[i] << endl;
  

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
