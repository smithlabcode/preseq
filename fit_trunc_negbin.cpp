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
  void estim_params(const vector<size_t> &vals);
  void estim_trunc_params_Newton(const vector<size_t> &vals_hist);

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

void NBD::estim_params(const vector<size_t> &vals){
  mu = accumulate(vals.begin(), vals.end(), 0.0)/vals.size(); 
  //mu= (1/n)sum(x), accumulate takes the sum of vals.begin
  
  // Now for the alpha
  const int max_value = *max_element(vals.begin(), vals.end()); //pointer to the max element of vals
  vector<int> vals_hist(static_cast<size_t>(max_value) + 1, 0.0);  //static_cast converts types
  for (size_t i = 0; i < vals.size(); ++i)
    ++vals_hist[static_cast<size_t>(vals[i])];  //vals_hist gives the histogram of the data ie how many x's = 1,2,...
  
  const double vals_count = vals.size(); //number of elements
  
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

static inline double
NBD_trunc_log_L_partial_w(const vector<size_t> &vals_hist,
                           const double w,
                           const double k){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0)); 
  const double mean = k*(1-w)/(w*(1-pow(w,k)));
  double log_L_partial_w = 
    vals_size*k/(w*(1-pow(w,k))) - vals_size*mean/(1-w);

  return(log_L_partial_w);
}

static inline double
NBD_trunc_log_L_partial_k(const vector<size_t> &vals_hist,
                          const double w,
                          const double k){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));

  double log_L_partial_k = 0;

  for(size_t i = 1; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += 1/(k+j);

    log_L_partial_k += vals_hist[i]*holding_val;
  }
  log_L_partial_k += vals_size*log(w)/(1-pow(w,k));
  return(log_L_partial_k);
}

static inline double 
NBD_trunc_log_L_2nd_partial_w(const vector<size_t> &vals_hist,
                               const double w,
                               const double k){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));
  const double denom = 1 - pow(w,k);
  const double mean = k*(1-w)/(w*denom);

  double log_L_2nd_partial_w = -vals_size*mean/((1-w)*(1-w)) -
    vals_size*k*(1 - (k+1)*pow(w,k))/(w*w*denom*denom);

  return(log_L_2nd_partial_w);
}

static inline double 
NBD_trunc_log_L_mixed_2nd_partial(const vector<size_t> &vals_hist,
                                 const double w,
                                 const double k){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));

  const double denom = 1 - pow(w,k);
  double log_L_mixed_2nd_partial = vals_size*(1 - (1-k*log(w))*pow(w,k))/
    (w*denom*denom);

  return(log_L_mixed_2nd_partial);
}

static inline double 
NBD_trunc_log_L_2nd_partial_k(const vector<size_t> &vals_hist,
                             const double w,
                             const double k){

  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));

  double log_L_2nd_partial_k = 0;
  const double denom = 1 - pow(w,k);

  for(size_t i = 1; i < vals_hist.size(); i++){
    double holding_val = 0;
    for(size_t j = 0; j < i; j++)
      holding_val += 1/((k+j)*(k+j));

    log_L_2nd_partial_k -= vals_hist[i]*holding_val;
  }
  log_L_2nd_partial_k += 
    vals_size*log(w)*log(w)*pow(w,k)/(denom*denom);

  return(log_L_2nd_partial_k);
}

static inline void
NBD_compute_inverse_Hessian(const vector<size_t> &vals_hist,
                            const double w,
                            const double k,
                            vector< vector<double> > &invHessian){

  const double log_L_2nd_partial_w = 
    NBD_trunc_log_L_2nd_partial_w(vals_hist,w,k);

  const double log_L_2nd_partial_k =
    NBD_trunc_log_L_2nd_partial_k(vals_hist,w,k);
  
  const double log_L_2nd_mixed_partial = 
    NBD_trunc_log_L_mixed_2nd_partial(vals_hist,w,k);

  const double denom = log_L_2nd_partial_w*log_L_2nd_partial_k
    - log_L_2nd_mixed_partial*log_L_2nd_mixed_partial;

  invHessian[0][0] = log_L_2nd_partial_k/denom;
  invHessian[1][1] = log_L_2nd_partial_w/denom;
  invHessian[1][0] = -log_L_2nd_mixed_partial/denom;
  invHessian[0][1] = invHessian[1][0];
}
   
   

void 
NBD::estim_trunc_params_Newton(const vector<size_t> &vals_hist){

  double sum_vals = 0.0;
  for(size_t i = 0; i < vals_hist.size(); i++)
    sum_vals += i*vals_hist[i];
  const double vals_size = 
    static_cast<double>(accumulate(vals_hist.begin(), vals_hist.end(), 0));
  const double mean = sum_vals/vals_size;
  double k_estim = 2;
  double w_estim = (k_estim+mean)/(mean+2*k_estim);

  double diff = numeric_limits<double>::max();
  double prev_log_L = numeric_limits<double>::max(); 
  double current_log_L = -1;

  vector< vector<double> > invHessian( 2, vector<double>(2,0));

  while(diff > tolerance){
    
    NBD_compute_inverse_Hessian(vals_hist,w_estim,k_estim,invHessian);
    double log_L_partial_w = 
      NBD_trunc_log_L_partial_w(vals_hist,w_estim,k_estim);
    double log_L_partial_k = 
      NBD_trunc_log_L_partial_k(vals_hist,w_estim,k_estim);

    w_estim += 0.01*(invHessian[0][0]*log_L_partial_w 
		     + invHessian[0][1]*log_L_partial_k);
    
    k_estim += 0.01*(invHessian[1][1]*log_L_partial_k +
		    invHessian[1][0]*log_L_partial_w);

    current_log_L = NBD_trunc_log_L(vals_hist,w_estim,k_estim);
    diff = fabs((current_log_L - prev_log_L)/max(current_log_L,prev_log_L));
    prev_log_L = current_log_L;
  
    double mu = k_estim*(1-w_estim)/w_estim;
    cerr << "mu = " << mu << ", " << "alpha = " << 1/k_estim << "\n";

  }
  mu = k_estim*(1-w_estim)/w_estim;
  alpha = 1/k_estim;
   

}

  

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;
    string resample_size_string;
    string mus_string;
    string alpha_string;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("fit_trunc_negbin", "",
			   "");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE); 

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

    NBD NBD_distro(1,2);
    NBD_distro.estim_trunc_params_Newton(vals_hist);
  
    

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




    
    

