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

class NBD {
public:
  NBD(const double m_, const double a_) : mu(m_), alpha(a_) {set_helpers();}

  double get_mu() const {return mu;}
  double get_alpha() const {return alpha;}
 
  void set_mu(const double m_) {mu = m_;}
  void set_alpha(const double a_) {alpha = a_;}
 
  double operator()(int val) const;
  void estim_params(const vector<int> &x);
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
  //TODO: check these!!!!!!!!!!!!!!!!  they are correct.
}

double
NBD::operator()(const int val) const {
  const double P = (gsl_sf_lngamma(val + n_helper) -
                    gsl_sf_lnfact(static_cast<size_t>(val))) +  //should this be val+1? no, it's calling the factorial function
    n_log_p_minus_lngamma_n_helper + val*log_q_helper;
  if (!finite(P)) return -40; //why return -40 if the log likelihood is -infty?
  return P;
}


static inline double
score_fun_first_term(const vector<int> &vals_hist,
                     const double mu, const double alpha) {
  double sum = 0;
  for (size_t i = 0; i < vals_hist.size(); ++i)
    if (vals_hist[i] > 0) {
      double inner_sum = 0;
      for (size_t j = 0; j < i; ++j)
        inner_sum += j/(1 + alpha*j);
      sum += vals_hist[i]*inner_sum;
    }
  return sum; //sum=sum_{i=1}^{max(x)-1} (#(x=i)) * (sum_{j=1}^{i-1} j/(1+alpha*j) )  but why? and where is mu?
}

static inline double
alpha_score_function(const vector<int> &vals_hist, const double mu,
                     const double alpha, const double vals_count) {  //vals_count= # x > 0
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_fun_first_term(vals_hist, mu, alpha)/vals_count +
          (log(one_plus_alpha_mu)/alpha - mu)/alpha);  //returns sum/(# x neq 0) + (log(1+alpha*mu)/alpha - mu)/alpha
}

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

void NBD::estim_params(const vector<int> &x){
  mu = accumulate(x.begin(), x.end(), 0.0)/x.size(); //mu= (1/n)sum(x), accumulate takes the sum of vals.begin
 
  // Now for the alpha
  const int max_value = *max_element(x.begin(), x.end()); //pointer to the max element of vals
  vector<int> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);  //static_cast converts types
  for (size_t i = 0; i < x.size(); ++i)
    ++vals_hist[static_cast<size_t>(x[i])];  //vals_hist gives the histogram of the data ie how many x's = 1,2,...
 
  const double vals_count = x.size(); //number of elements
 
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
 
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
 
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_count);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  //bisection, but what happened to the terms involving the gamma func? See Zhang et al. top of page 7
  set_helpers();
}

void
generate_mix_NBD(const vector<double> &mixing,
                 const vector<NBD> &distros,
                 vector<size_t> &sample){

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
        const double n = 1/distros[j].get_alpha();
        const double p = n/(n + distros[j].get_mu());
        sample[i] = gsl_ran_negative_binomial(rng, p, n);
        break;
      }
    }
  }
     
}

//generate obs_probs according to a dirichlet dist with gamma priors
void
generate_obs_probs(const vector<double> &thetas,
                   const vector<double> &alphas,
                   const vector<double> &mus,
                   vector<double> &observation_probabilities){

  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());
  
  //generate gamma priors
  double prior_weights[observation_probabilities.size()];
  for(size_t i = 0; i < observation_probabilities.size(); i++){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < thetas.size(); i++){
      test_val += thetas[j];
      if(u <= test_val){
        const double a = 1/alphas[j];
        const double b = mus[j]*alphas[j];
        prior_weights[i] = gsl_ran_gamma(rng, a, b);
        break;
      }
    }
  }
        
  double holding_probs[observation_probabilities.size()];
  gsl_ran_dirichlet(rng, observation_probabilities.size(),
                    prior_weights, holding_probs);
  for(size_t i = 0; i < observation_probabilities.size(); i++)
    observation_probabilities[i] = holding_probs[i];
}    

//generate reads from probabilities
void
generate_reads(const vector<double> &observation_probabilities,
               const size_t number_reads,
               vector<size_t> &read_counts){
  
  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());
  
  for(size_t i = 0; i < number_reads; i++){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < read_counts.size(); j++){
      test_val += observation_probabilities[j];
      if(u <= test_val){
        read_counts[j]++;
        break;
      }
    }
  }
}

double
generate_reads_until_epsilon_reads(const vector<double> &observation_probabilities,
                                            const double epsilon,
                                            vector<size_t> &read_counts){
  double hitting_time = 0.0;
  
  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());
  
  double prop_obs_reads =0.0;
  double transcript_weight = 1/static_cast<double>(read_counts.size());
  
  while(prop_obs_reads < 1- epsilon){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < read_counts.size(); j++){
      test_val += observation_probabilities[j];
      if(u <= test_val){
        if(read_counts[j] == 0)
          prop_obs_reads += transcript_weight;
        read_counts[j]++;
        hitting_time++;
        break;
      }
    }
  }
  return(hitting_time);
}

double
generate_reads_until_all_epsilon_mostlikely(const vector<double> &observation_probabilities,
                                            const double epsilon,
                                            vector<size_t> &read_counts){
  
  size_t hitting_count = floor((1-epsilon)*read_counts.size());
  cerr <<"hitting count = " << hitting_count << "\n";
  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());
  size_t count = 0;
  double hitting_time = 0;
  
  while(count < hitting_count){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < read_counts.size(); j++){
      test_val += observation_probabilities[j];
      if(u <= test_val){
        if(read_counts[j] == 0 && j >= hitting_count)
          count++;
        
        read_counts[j]++;
        hitting_time++;
        break;
      }
    }
  }
  return(hitting_time);
}

double
generate_reads_until_delta_all_epsilon_mostlikely(const vector<double> &observation_probabilities,
                                                  const double epsilon,
                                                  const double delta,
                                                  vector<size_t> &read_counts){
  
  size_t hitting_count = floor((1-delta)*(1-epsilon)*read_counts.size());
  const gsl_rng_type *T;
  gsl_rng *rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL) + getpid());
  size_t count = 0;
  double hitting_time = 0;
  
  while(count < hitting_count){
    double u = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < read_counts.size(); j++){
      test_val += observation_probabilities[j];
      if(u <= test_val){
        if(read_counts[j] == 0 && j >= hitting_count)
          count++;
        
        read_counts[j]++;
        hitting_time++;
        break;
      }
    }
  }
  return(hitting_time);
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    size_t number_transcripts = 1000000;
    size_t number_mixtures = 1;
    string alphas_string;
    string mus_string;
    string thetas_string;
    size_t number_reads = 10000;
    double epsilon = 0.05;
    double delta = 0.05;
    size_t bootstraps = 1000;
   
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("negbin_generator", "",
                           "");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false , VERBOSE);
    opt_parse.add_opt("number_transcripts", 't', "total number transcripts",
                      false , number_transcripts);
    opt_parse.add_opt("mixing", 'x', "mixing parameters, separated by commas",
                      true, thetas_string);
    opt_parse.add_opt("number_mixtures", 'm', "number of mixtures",
                      false, number_mixtures);
    opt_parse.add_opt("alphas",'a',"alphas, seperated by commas",
                      true, alphas_string);
    opt_parse.add_opt("mus",'u',"mus, seperated by commas",
                      true, mus_string);
    opt_parse.add_opt("number_reads",'n',"Initial number of reads to simulate",
                      false, number_reads);
    opt_parse.add_opt("epsilon",'e',"Epsilon",
                      false, epsilon);
    opt_parse.add_opt("delta",'d',"Delta",
                      false, delta);
    opt_parse.add_opt("bootstraps",'b',"number of bootstraps",
                      false, bootstraps);
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

    vector<double> thetas;
    vector<double> alphas;
    vector<double> mus;

    vector<string> alphas_parts = rmap::split(alphas_string, ",");    
    vector<string> mus_parts  = rmap::split(mus_string, ",");                                                     
    vector<string> thetas_parts = rmap::split(thetas_string, ",");



    for(size_t i = 0; i < number_mixtures; i++){
      alphas.push_back(atof(alphas_parts[i].c_str()));
      mus.push_back(atof(mus_parts[i].c_str()));
      thetas.push_back(atof(thetas_parts[i].c_str()));
    }
    double mix_sum = accumulate(thetas.begin(), thetas.end(), 0.0);
    if( mix_sum!= 1){
      cerr << "Mixing parameters must sum to 1! \n";
      cerr << "Mixing sums to " <<
        mix_sum << ", need to reweight";
    }
    for(size_t i = 0; i < thetas.size(); i++)
      thetas[i] = thetas[i]/mix_sum;


    vector<double> observation_probabilities(number_transcripts, 0.0);

    generate_obs_probs(thetas, alphas, mus, observation_probabilities);
    sort(observation_probabilities.begin(), observation_probabilities.end());
    

    vector<double> hitting_times;
    for(size_t i = 0; i < bootstraps; i++){
      vector<size_t> read_counts(number_transcripts, 0);

      double hitting_time = generate_reads_until_delta_all_epsilon_mostlikely(observation_probabilities,
                                                                      epsilon, delta, read_counts);
      hitting_times.push_back(hitting_time);
    }
    
    ostream* out = (outfile.empty()) ?
    &std::cout : new std::ofstream(outfile.c_str());
    for(size_t i = 0; i < hitting_times.size(); i++){
      *out << hitting_times[i] << endl;
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
     
