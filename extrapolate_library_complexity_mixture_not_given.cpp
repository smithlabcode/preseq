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


static void
sample_genome_count_reads(const gsl_rng *rng,
                          const vector<size_t> &genome_lengths,
                          const vector<double> mixing,
                          const size_t step_size,
                          vector<bool> &the_sample,
                          size_t &total_reads,
                          size_t &unique_reads){

 
  for(size_t i = 0; i < step_size; i++){

    size_t holding_val = 0;
    const double u = gsl_rng_uniform(rng);
    const double v = gsl_rng_uniform(rng);
    double test_val = 0.0;
    for(size_t j = 0; j < mixing.size(); j++){
      test_val += mixing[j];

      if(u < test_val){
        holding_val += 
          static_cast<size_t>(floor(v*genome_lengths[j]));
        break;
      }

      holding_val += genome_lengths[j];
    }

    assert(holding_val <= the_sample.size());

    if(the_sample[holding_val] == 0){
      ++unique_reads;
      the_sample[holding_val] = 1;
    }
    ++total_reads;
  }
}

static inline size_t
get_max(const size_t x1, const size_t x2){
  size_t larger = (x1 > x2) ? x1 : x2;
  return(larger);
}

static inline size_t
calculate_zeros(const double lambda,
                const double mixing,
                const size_t values_size){
  return(ceil(mixing*values_size*exp(-lambda)/(1-exp(-lambda))));
}

static void
calculate_genome_lengths(const vector<double> &lambdas,
                         const size_t truncation_max,
                         const size_t values_size,
                         vector<size_t> &genome_lengths,
                         vector<double> &mixing){
  double mean  = 0;
  for(size_t i = 0; i < lambdas.size(); i++)
    mean += mixing[i]*lambdas[i]/(1-exp(-lambdas[i]));
  
  vector<size_t> observed(lambdas.size(), 0);
  for(size_t i = 0; i < lambdas.size(); i++)
    observed[i] = ceil(mixing[i]*values_size);

  for(size_t i = 0; i < genome_lengths.size(); i++){
    genome_lengths[i] = observed[i];
    genome_lengths[i] += calculate_zeros(lambdas[i], mixing[i],
                                         values_size);
    if(genome_lengths[i] > truncation_max)
      genome_lengths[i] = truncation_max;
  }
  /* This is the previous method of calculating genome length

  for(size_t i = 0; i < genome_lengths.size(); i++){
    genome_lengths[i] = round(observed[i]/lambdas[i]);
    genome_lengths[i] = get_max(observed[i], genome_lengths[i]);
  }

  G = N/lambda
  */

  /* mixing should be redefined to relative mixing, ie
     number of hits in group j/ total hits 
     = relative mean/total mean  */
  for(size_t i = 0; i < mixing.size(); i++)
    mixing[i] = mixing[i]*lambdas[i]/((1-exp(-lambdas[i]))*mean);
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    size_t step_size = 1000000;
    size_t resample_size = 10000000;
    size_t truncation_max = 1e9;

    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("extrapolate_library_complexity", "",
			   "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("step_size", 'p', "step size between resamples",
                      false, step_size);
    opt_parse.add_opt("resample_size",'r',"total sample size of resampled data",
                      false, resample_size);
    opt_parse.add_opt("truncation_max", 't', "max length of individual genomes",
                      false, truncation_max);

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
    vector<double> mixing;

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
      lambdas.push_back(holding_val);
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
      mixing.push_back(holding_val);
      indx++;
    }
    in.close();
 
    vector<size_t> genome_lengths(number_mixtures, 0);
    calculate_genome_lengths(lambdas, truncation_max, values_size, 
                             genome_lengths, mixing);

    size_t unique_reads = 0;
    size_t total_reads = 0;
    size_t total_genome_length = accumulate(genome_lengths.begin(),
                                            genome_lengths.end(), 0);
    vector<bool> resampled_sample(total_genome_length, 0);

 
    const gsl_rng_type *T;
    gsl_rng *rng;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, time(NULL) + getpid());

    
 
    ostream* out = (outfile.empty()) ? 
      &std::cout : new std::ofstream(outfile.c_str());
    *out << 0 << "\t" << 0 << endl;
    for(size_t i = 0; i < resample_size; i += step_size){
      sample_genome_count_reads(rng, genome_lengths, mixing, 
                                step_size, resampled_sample, 
                                total_reads, unique_reads);
      *out << total_reads << "\t" << unique_reads << endl;
    }
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
