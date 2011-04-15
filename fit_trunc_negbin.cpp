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

static void
sample_vals(const gsl_rng *rng,
            const size_t resample_size,
            const vector<size_t> &orig,
            vector<size_t> &new_vals){
  vector<size_t> new_vals_indx(resample_size, 0);
  gsl_ran_choose(rng, (size_t *)&new_vals_indx.front(), resample_size,
                 (size_t *)&orig.front(), orig.size(), sizeof(size_t));

  size_t count = 1;
  for(size_t i = 1; i < new_vals_indx.size(); i++){
    if(new_vals_indx[i] == new_vals_indx[i-1])
      count++;
    else{
      new_vals.push_back(count);
      count = 1;
    }
  }
  new_vals.push_back(count);
}
                
    


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
  void estim_trunc_params_Newton(const vector<size_t> &vals_hist);
  double trunc_log_L(const vector<size_t> &vals_hist);
  double trunc_log_pdf(size_t val);
 

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
score_fun_first_term(const vector<size_t> &vals_hist, 
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
alpha_score_function(const vector<size_t> &vals_hist, const double mu, 
		     const double alpha, const double vals_count) {  //vals_count= # x > 0
  const double one_plus_alpha_mu = 1 + alpha*mu;
  return (score_fun_first_term(vals_hist, mu, alpha)/vals_count + 
	  (log(one_plus_alpha_mu)/alpha - mu)/alpha);  //returns sum/(# x neq 0) + (log(1+alpha*mu)/alpha - mu)/alpha
}

static inline double
movement(const double a, const double b) {
  return fabs((a - b)/max(a, b)); //delta
}

static double
compute_mean(const vector<size_t> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(),
				   vals_hist.end(), 0));
  double mean = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}
    

void 
NBD::estim_params(const vector<size_t> &vals_hist){
  mu = compute_mean(vals_hist);
  //mu= (1/n)sum(x), accumulate takes the sum of vals.begin
  
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), 
				   vals_hist.end(), 0));
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(vals_hist, mu, a_mid, vals_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  //bisection, but what happened to the terms involving the gamma func? See Zhang et al. top of page 7
  set_helpers();
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


static double
compute_mean(const vector<double> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(),
				   vals_hist.end(), 0.0));
  double mean = 0.0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    mean += i*vals_hist[i]/vals_size;
  }
  return(mean);
}
    

void 
NBD::estim_params(const vector<size_t> &vals_hist,
		  const vector<double> &probs){
  vector<double> pseudo_hist(vals_hist.size(), 0.0);
  for(size_t i = 0; i < vals_hist.size(); i++){
    pseudo_hist[i] = vals_hist[i]*probs[i];
  }
  mu = compute_mean(pseudo_hist);
  
  const double pseudo_size = 
    accumulate(pseudo_hist.begin(), pseudo_hist.end(), 0.0);
  
  double a_low = min_allowed_alpha;
  double a_high = max_allowed_alpha;
  double a_mid = max_allowed_alpha;
  double mid_val;
  
  double diff = numeric_limits<double>::max();
  double prev_val = numeric_limits<double>::max();
  
  while (diff > tolerance && movement(a_high, a_low) > tolerance) {
    a_mid = (a_low + a_high)/2;
    mid_val = alpha_score_function(pseudo_hist, mu, a_mid, pseudo_size);
    if (mid_val < 0) a_high = a_mid;
    else a_low = a_mid;
    diff = fabs((prev_val - mid_val)/prev_val);
    prev_val = mid_val;
  }
  alpha = a_mid;  //bisection, but what happened to the terms involving the gamma func? See Zhang et al. top of page 7
  set_helpers();
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

double
NBD::trunc_log_L(const vector<size_t> &vals_hist){
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));
  double log_L = 0;
  const double k = 1/alpha;
  const double w = k/(mu+k);
  for(size_t i = 1; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += log(k+j);

    log_L += vals_hist[i]*gsl_sf_lngamma(i) 
                        + vals_hist[i]*holding_val;

  }
 
  log_L += vals_size*k*log(w) - vals_size*log(1-pow(w,k));
    
  return(log_L);
}

/* k = 1/alpha, w = k/(mu+k), following sampford(Biometrika 1955) */
static inline double
NBD_trunc_log_L(const vector<size_t> &vals_hist,
                const double w,
                const double k){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));
  double log_L = 0;
  for(size_t i = 1; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += log(k+j);

    log_L += vals_hist[i]*gsl_sf_lngamma(i) 
                        + vals_hist[i]*holding_val;

  }
 
  log_L += vals_size*k*log(w) - vals_size*log(1-pow(w,k));
    
  return(log_L);
}

double 
NBD::trunc_log_pdf(const size_t val){
  const double k = 1/alpha;
  const double w = k/(mu+k);
  double holding_val = 0.0;
  if(val > 0){
    for(size_t j =0; j < val; ++j){
      holding_val += log(k+j);
    }
    holding_val += gsl_sf_lngamma(val);
  }
  return(holding_val + k*log(w) - log(1-pow(w,k)) + val*log(1-w));
}

static void
nonzero_expectation_step(const vector<size_t> &vals_hist,
                         vector<NBD> &distros,
                         const vector<double> &mixing,
                         vector< vector<double> > &probs){
  for(size_t i = 0; i < vals_hist.size(); i++){
    vector<double> log_denom_vec;
    
    for(size_t j = 0; j < distros.size(); j++){
      log_denom_vec.push_back(log(mixing[j]) +
                              distros[j].trunc_log_pdf(i));
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

  }
}


static double
NBD_add_on_expected_zeros(const NBD &distro,
			  const double pseudo_size){
  const double alpha = distro.get_alpha();
  const double mu = distro.get_mu();
  const double prob_zero = pow(1+alpha*mu, -1/alpha);
  const double expected_zeros = pseudo_size*(prob_zero/(1-prob_zero));
  return(expected_zeros);
}


static void
nonzero_maximization_step(const vector<size_t> &vals_hist,
                          const vector< vector<double> > &probs,
                          vector<NBD> &distros){

  for(size_t i = 0; i < distros.size(); i++){
    distros[i].estim_params(vals_hist, probs[i]);
  }
}

static double 
trunc_log_L_hist(const vector<size_t> &vals_hist,
                 vector<NBD> &distros,
                 const vector<double> &mixing){
  double logL = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++){
    double inner_sum = 0.0;
    vector<double> workspace_vec(distros.size(), 0.0);
    for(size_t j = 0; j < distros.size(); j++)
      workspace_vec[j] = log(mixing[j]) + 
                         distros[j].trunc_log_pdf(i);
    
    inner_sum = log_sum_log_vec(workspace_vec, workspace_vec.size());
    logL += vals_hist[i]*inner_sum;
  }
  return(logL);
}

double 
NBD_EM_add_zeros_mixture(const double &tol,
			 const size_t max_iter,
			 vector<size_t> &vals_hist,
			 vector<NBD> &distros,
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

    nonzero_expectation_step(vals_hist, distros, 
                                     mixing, probs);
 
    double expected_zeros = 0.0;
    for(size_t j = 0; j < distros.size(); j++){
      double psuedo_size = 0.0;
      for(size_t l = 1; l < vals_hist.size(); l++){
	psuedo_size += vals_hist[l]*probs[j][l];
      }
      expected_zeros += NBD_add_on_expected_zeros(distros[j], psuedo_size);
    }
    vals_hist[0] = round(expected_zeros);

    nonzero_maximization_step(vals_hist, probs, distros);
    vals_hist[0] = 0;
    calculate_mixing(vals_hist, probs, mixing);

    score = trunc_log_L_hist(vals_hist, distros, mixing);
    error = fabs((score - prev_score)/score);
    if(error < tol){
      break;
    }
    prev_score = score;

  }

  return(trunc_log_L_hist(vals_hist, distros, mixing));
}

static inline double
get_max(const double x1, const double x2){
  const double larger = (x1 > x2) ? x1 : x2;
  return(larger);
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;
    size_t number_mixtures = 20;
    size_t max_iter = 100000;
    double tolerance = 1e-20;
    size_t starting_number_mixtures = 1;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("fit_trunc_negbin", "",
			   "");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("number_mixtures", 'm', 
                      "number of mixtures to take",
                      false, number_mixtures);
    opt_parse.add_opt("tolerance", 't', "Numerical tolerance",
                     false, tolerance);
    opt_parse.add_opt("starting_number_mixtures",'s',"Starting # mixtures",
		      false, starting_number_mixtures);

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
    const string input_file_name = leftover_args.front();

    /*****************************************************************/

   

    vector<size_t> vals;
    std::ifstream in(input_file_name.c_str());
    if (!in) 
      throw BEDFileException("cannot open input file " + input_file_name);
    static const size_t buffer_size = 10000; // Magic!
    in.peek();
    while (!in.eof()){
      int holding_val = 0;
      char buffer[buffer_size];
      in.getline(buffer, buffer_size);
      holding_val = atoi(buffer);
      if(holding_val > 0)
        vals.push_back(holding_val);
    }
    in.close();

    const size_t max_value = *std::max_element(vals.begin(),vals.end());
    vector<size_t> vals_hist(max_value + 1, 0.0);
    for (size_t i = 0; i < vals.size(); ++i){
      ++vals_hist[static_cast<size_t>(vals[i])];
    }


      vector<NBD> distros;
      vector<double> mixing(number_mixtures, 
                          1/static_cast<double>(number_mixtures));
    
      sort(vals.begin(), vals.end());
      double step_size = floor(static_cast<double>(vals.size())
                                 /static_cast<double>(number_mixtures));
      for(size_t j = 0; j < number_mixtures; j++){
        double mu_hat = accumulate(vals.begin()+j*step_size,
                              vals.begin()+(j+1)*step_size, 
                              0.0)/step_size;

        if(j != 0){
          NBD last_distro = distros.back();
          mu_hat = get_max(mu_hat, last_distro.get_mu()) 
                         + 0.01;
        }
        else
          mu_hat = mu_hat+0.01;

 
        NBD holding_distro(mu_hat,1);
        distros.push_back(holding_distro);
      }
      double score = 0.0;  
      score  = NBD_EM_add_zeros_mixture(tolerance, max_iter, vals_hist,
                                   distros, mixing);
     ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    *out << vals.size() << endl;
    *out << number_mixtures << endl;
    *out << score - (2*number_mixtures - 1) << endl;
    for(size_t j = 0; j < distros.size(); j++){
      *out << distros[j].get_mu() << "\t" << distros[j].get_alpha() << endl;
    } 

    for(size_t j = 0; j < mixing.size(); j++){
      *out << mixing[j] << endl;
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




    
    

