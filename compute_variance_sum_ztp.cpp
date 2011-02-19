#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>


#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])


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

static void
compute_log_normalizing_constant(const double lambda,
                             double &log_normalizing_constant){
    log_normalizing_constant -= log(exp(lambda) - 1);
}


static void
compute_current_lambda_vec(const double current_lambda,
                           const size_t upper_lim,
                           vector<double> &normalized_lambdas){

  normalized_lambdas[0] = current_lambda;
  for(size_t i = 1; i < upper_lim; i++)
    normalized_lambdas[i] = normalized_lambdas[i-1]*current_lambda/(i+1);

}
  
/*correct the convolution */
static void
set_initial_Qs_to_zero(const size_t iter,
                       vector<double> &current_Q){
  for(size_t i = 0; i < iter; i++)
    current_Q[i] = 0;
}

static void
set_back_Qs_to_zero(const size_t upper_lim,
                    vector<double> &current_Q){
  for(size_t i = upper_lim; i < current_Q.size(); i++)
	current_Q[i] = 0.0;
}

void
compute_current_Q_by_fft(const size_t upper_lim2,
	                 vector<double> &prev_Q,
                         vector<double> &normalized_lambdas,
                         vector<double> &current_Q){
  /*I like to use vectors, but gsl_FFt does not */
  double lambdas_array[2*upper_lim2];
  double prev_Q_array[2*upper_lim2];
  double current_Q_array[2*upper_lim2];
  for(size_t i = 0; i < upper_lim2; i++){
    REAL(lambdas_array, i) = normalized_lambdas[i];
    IMAG(lambdas_array, i) = 0.0;
    REAL(prev_Q_array, i) = prev_Q[i];
    IMAG(prev_Q_array, i) = 0.0;
    REAL(current_Q_array, i) = 0.0;
    IMAG(current_Q_array, i) = 0.0;
  }

  gsl_complex_packed_array gsl_lambdas = lambdas_array;
  gsl_complex_packed_array gsl_prev_Q = prev_Q_array;
  gsl_complex_packed_array gsl_current_Q = current_Q_array;
  /*fft by gsl */
  gsl_fft_complex_radix2_forward(gsl_lambdas, 1, upper_lim2); 
  gsl_fft_complex_radix2_forward(gsl_prev_Q, 1, upper_lim2);
  /*multiply fft(Q) and fft(lambdas) to get fft(convolution) */
  for(size_t i = 0; i < upper_lim2; i++){
    REAL(gsl_current_Q, i) = REAL(gsl_prev_Q, i)*REAL(gsl_lambdas, i) -
      IMAG(gsl_prev_Q, i)*IMAG(gsl_lambdas, i);
    IMAG(gsl_current_Q, i) = REAL(gsl_prev_Q, i)*IMAG(gsl_lambdas, i) +
      REAL(gsl_lambdas, i)*IMAG(gsl_prev_Q, i);
  }
  /*invert fft to set current_Q = convolution*/
  gsl_fft_complex_radix2_inverse(gsl_current_Q, 1, upper_lim2);
  for(size_t i = 0; i < upper_lim2; i++){
    current_Q[i] = REAL(gsl_current_Q, i);
  }
  
}

/*need to invert rows and columns */ 
static void
compute_cdf_sum_ztp(const vector< vector<double> > &cdf_sum_ztp,
                    vector< vector<double> > &cdf_inverse_sum_ztp){

  for(size_t i = 0; i < cdf_inverse_sum_ztp.size(); i++){
    for(size_t j = 0; j < cdf_inverse_sum_ztp[i].size(); j++){
      cdf_inverse_sum_ztp[i][j] = 1 - exp(cdf_sum_ztp[j][i]);
    }
  }
}
  








int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;
    size_t upper_lim = 10000000;
    size_t step_size = 1000000;

    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("compute_variance_sum_stp", "",
			   "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("upper limit",'u',"upper limit of n",
                      false, upper_lim);
    opt_parse.add_opt("step_size", 's', "step size to take probabilities",
                      false, step_size);


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

    vector<double> lambdas;

    std::ifstream in(input_file_name.c_str());
    if (!in) 
      throw BEDFileException("cannot open input file " 
                               + input_file_name);

    static const size_t buffer_size = 10000; // Magic!
    while (!in.eof()){
      double holding_val = 0;
      char buffer[buffer_size];
      in.getline(buffer, buffer_size);
      holding_val = atof(buffer);
      if(!(holding_val >= 0)){
        cerr << "Invalid Input.\n";
      }
      assert(holding_val >= 0);
      lambdas.push_back(holding_val);
    }
    in.close();
    if(lambdas.back() == 0)
      lambdas.pop_back();

    size_t upper_lim2 = 1;
    while((static_cast<double>(upper_lim2)/static_cast<double>(upper_lim))
          < 1){
      upper_lim2 = upper_lim2*2;
    }
    vector<double> prev_Q(upper_lim2, 0.0);
    prev_Q[0] = 1;

    size_t row_size = upper_lim/step_size;
    vector< vector<double> > cdf_sum_ztp(lambdas.size(), 
			       vector<double>(row_size ,0.0));

    double log_normalizing_constant = 0.0;
    size_t column_count = 0;
    
    for(size_t i = 0; i < lambdas.size(); i++){
      vector<double> normalized_lambda_vec(upper_lim2, 0.0);
      vector<double> current_Q(upper_lim2, 0.0);

      compute_current_lambda_vec(lambdas[i], upper_lim, normalized_lambda_vec);
      compute_current_Q_by_fft(upper_lim2, prev_Q,
                               normalized_lambda_vec, current_Q);
      set_initial_Qs_to_zero(column_count, current_Q);

      set_back_Qs_to_zero(upper_lim, current_Q);
      
      compute_log_normalizing_constant(lambdas[i], log_normalizing_constant);
      for(size_t j = 0; j < row_size; j++){
        cdf_sum_ztp[column_count][j] = log(accumulate(current_Q.begin(),
	  current_Q.begin()+(j+1)*step_size -1, 0.0));
        cdf_sum_ztp[column_count][j] += log_normalizing_constant;
        
        
      }
      column_count++;

      prev_Q.swap(current_Q);
    }
    vector< vector<double> > cdf_inverse_sum_ztp(row_size, 
                                         vector<double>(lambdas.size(), 0.0));
    compute_cdf_sum_ztp(cdf_sum_ztp, cdf_inverse_sum_ztp);

    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    for(size_t i = 0; i < cdf_inverse_sum_ztp.size(); i++){
      *out << i+1 << "\t";
      for(size_t j = 0; j < cdf_inverse_sum_ztp[i].size(); j++)
        *out << cdf_inverse_sum_ztp[i][j] << "\t";
      
      *out << endl;
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


    
      
