/*    etcomplexity
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *                       Timothy Daley
 *
 *    Authors: Andrew D. Smith and Timothy Daley
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
#include "MappedRead.hpp"

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>

#include <gsl/gsl_vector_double.h>  
#include <gsl/gsl_matrix_double.h>  
#include <gsl/gsl_linalg.h>  

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;

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
      // abort if the sum is infinte //
      assert(finite(sum)); 
#endif
    }
  }
  return(max_val + log(sum));
}


static inline double
weight_exponential(const double dist, double decay_factor) {
  return std::pow(0.5, decay_factor*dist);
}


static void
smooth_histogram(const size_t bandwidth,
		 const double decay_factor, vector<double> &hist) {
  vector<double> updated_hist(hist);
  for (size_t i = 0; i < hist.size(); ++i) {
    double total_prob = 0.0, total_weight = 0.0;
    for (size_t j = ((i >= bandwidth/2) ? i - bandwidth/2 : 0);
	 j < std::min(i + bandwidth/2, hist.size()); ++j) {
      const double dist = std::abs(int(i) - int(j));
      const double w = weight_exponential(dist, decay_factor);
      total_prob += w*hist[j];
      total_weight += w;
    }
    updated_hist[i] = total_prob/total_weight;
  }
  updated_hist.swap(hist);
}


static void
solve_linear_system(const vector<vector<double> > &U, 
		    const vector<double> &v, vector<double> &b) {
  
  // compute Ub = v
  
  gsl_vector *v_gsl = gsl_vector_alloc(v.size());
  copy(v.begin(), v.end(), v_gsl->data);
  
  gsl_vector *b_gsl = gsl_vector_alloc(U.front().size());
  gsl_vector_set_zero(b_gsl);
  
  gsl_matrix *U_gsl = gsl_matrix_alloc(U.size(), U.front().size());
  for (size_t i = 0; i < U.size(); ++i)
    for (size_t j = 0; j < U[i].size(); ++j)
      gsl_matrix_set(U_gsl, i, j, U[i][j]);
  
  gsl_linalg_HH_solve(U_gsl, v_gsl, b_gsl);
  gsl_matrix_free(U_gsl);
  gsl_vector_free(v_gsl);
  
  copy(b_gsl->data, b_gsl->data + b_gsl->size, back_inserter(b));
  gsl_vector_free(b_gsl);
}



static void
compute_denom_ceoffs(const vector<double> &coeffs, const size_t denom_terms,
		     vector<double> &denom_coeffs) {
  const size_t m = denom_terms - 1;
  const size_t n = denom_terms;
  
  vector<double> v(denom_terms, 0.0);
  for (size_t j = 0; j < denom_terms; j++)
    v[j] = -coeffs[j + m + 1];
  
  //   for(j = m + 1; j < m + n + 1; j++){  
  //     gsl_vector_set(v, j - m - 1, - gsl_vector_get(c, j));  
  //     for(k = 1; k <= (j > n ? n : j); k++){  
  //       gsl_matrix_set(U, j - m - 1, k - 1, gsl_vector_get(c, j - k));  
  //     }  
  //   }  
  
  vector<vector<double> > U(n, vector<double>(n));
  for (size_t j = m + 1; j < m + n + 1; j++)
    for (size_t k = 1; k <= (j > n ? n : j); k++)
      U[j - m - 1][k - 1] = coeffs[j - k];
  
  solve_linear_system(U, v, denom_coeffs);
}


static void
test_coefficients(const vector<double> &c, const vector<double> &a, 
		  const vector<double> &b) {

  static const double COEFF_TEST_TOLERANCE = 1e-10;
  
  for (size_t i = 0; i < a.size(); ++i) {
    double sum = c[i];
    for (size_t j = 0; j < i; ++j)
      sum += b[j]*c[i - j - 1];
    assert(a[i] == sum);
  }

  for (size_t i = 0; i < b.size(); ++i) {
    const size_t offset = b.size() + i;
    double sum = c[offset];
    for (size_t j = 0; j < b.size(); ++j)
      sum += c[offset - 1 - j]*b[j];
    assert(sum < COEFF_TEST_TOLERANCE);
  }
}


static void
compute_pade_coeffs(const vector<double> &coeffs, const size_t num_terms,
		    vector<double> &num_coeffs, vector<double> &denom_coeffs) {
  
  const size_t denom_terms = coeffs.size() - num_terms;
  
  assert(num_terms + denom_terms == coeffs.size());
  
  denom_coeffs.clear();
  compute_denom_ceoffs(coeffs, denom_terms, denom_coeffs);
  
  num_coeffs.clear();
  num_coeffs.resize(num_terms, 0.0);
  for (size_t i = 0; i < num_coeffs.size(); ++i) {
    num_coeffs[i] = coeffs[i];
    for (size_t j = 0; j < std::min(i, denom_terms); ++j)
      num_coeffs[i] += denom_coeffs[j]*coeffs[i - j - 1];
  }
  
  test_coefficients(coeffs, num_coeffs, denom_coeffs);

}


static double
compute_pade_approx_numerator(const double t,
			      const vector<double> &numerator) {
  //   unsigned int j;  
  //   for(j = 0; j < a->size; j++)  
  //     num += gsl_vector_get(a, j) * pow(x, j);  
  double total = 0.0, prev_total = 0.0;
  for (size_t j = 0; j < numerator.size(); ++j) {
    total += numerator[j]*pow(t, j);
    // CHECK FOR CANCELLATION!!!!!
    assert(total != prev_total);
    prev_total = total;
  }
  return total;
}


static double
compute_pade_approx_denominator(const double t,
				const vector<double> &denominator) {
  //   for(j = 0; j <= b->size; j++)  
  //     den += (j == 0 ? 1.0 : gsl_vector_get(b, j - 1)) * pow(x, j);  
  double total = 1.0, prev_total = 1.0;
  for(size_t j = 0; j < denominator.size(); ++j) {
    total += denominator[j]*pow(t, j + 1);
    // CHECK FOR CANCELLATION!!!!!
    assert(total != prev_total);
    prev_total = total;
  }
  return total;
}


static bool 
check_sorted(const vector<MappedRead> &mr) {
  for (size_t i = 1; i < mr.size(); ++i)
    if (mr[i].r < mr[i-1].r)
      return false;
  return true;
}


static void
get_read_multiplicities(const vector<MappedRead> &reads, 
			vector<size_t> &values) {
  values.push_back(1);
  for (size_t i = 1; i < reads.size(); ++i) {
    if (reads[i].r == reads[i-1].r)
      values.back()++;
    else values.push_back(1);
  }
}


int
main(int argc, const char **argv){

  try {
    
    static const size_t MAGIC_10 = 10ul;
    
    /* FILES */
    string outfile;
    bool VERBOSE = false;
    
    size_t smoothing_bandwidth = 4;
    double smoothing_decay_factor = 8.0;
    
    double end_time = 10.0;
    double start_time = 2.0;
    double time_increment = 0.01;
    size_t max_terms = 10;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(argv[0], "", "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("lim", 'l', "extrapolation end time", false, end_time);
    opt_parse.add_opt("lower", 'L', "extrapolation start time", false, 
		      start_time);
    opt_parse.add_opt("inc", 'i', "time increment", false, time_increment);
    opt_parse.add_opt("terms",'T',"maximum terms to use", false, max_terms);
    
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
    
    if (VERBOSE)
      cerr << "[STATUS] loading reads ";
    
    std::ifstream in(input_file_name.c_str());
    MappedRead mr;
    vector<MappedRead> reads;
    while (in >> mr)
      reads.push_back(mr);
    if (VERBOSE)
      cerr << "[DONE]" << endl;
    if (reads.empty())
      throw RMAPException("could not get reads from: " + input_file_name);
    if (!check_sorted(reads))
      throw RMAPException("reads not sorted: " + input_file_name);
    
    vector<size_t> values;
    get_read_multiplicities(reads, values);
    
    // size_t tmp = 0ul;
    // while (in >> tmp)
    //   values.push_back(tmp);
    
    const size_t total_reads = accumulate(values.begin(), values.end(), 0ul); 
    // reads.size();
    if (VERBOSE)
      cerr << "TOTAL READS (N): " << total_reads << endl;
    
    vector<MappedRead>().swap(reads);
    
    const size_t max_multiplicity = 
      *std::max_element(values.begin(), values.end());
    const size_t distinct_reads = values.size();
    if (VERBOSE)
      cerr << "DISTINCT READS: " << distinct_reads << endl
	   << "PROB DISTINCT: " 
	   << static_cast<double>(distinct_reads)/total_reads << endl
	   << "MAX MULT: " << max_multiplicity << endl;
    
    vector<double> histogram(max_multiplicity, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++histogram[values[i] - 1];
    
    if (VERBOSE) {
      cerr << "MULTIPLICITIES WITH FREQ > 0: " 
	   << count_if(histogram.begin(), histogram.end(),
		       std::bind2nd(std::greater<double>(), 0.0)) << endl;
      cerr << "FIRST " << MAGIC_10 << " HIST VALUES:" << endl;
      for (size_t i = 0; i < MAGIC_10; ++i)
	cerr << i + 1 << "\t" << histogram[i] << endl; 
      cerr << endl;
    }

    smooth_histogram(smoothing_bandwidth, smoothing_decay_factor, histogram);
    
    vector<double> coeffs;
    copy(histogram.begin(), histogram.begin() + max_terms, 
	 back_inserter(coeffs));

    size_t non_zero_terms = 0;
    while (non_zero_terms < coeffs.size() && 
	   coeffs[non_zero_terms] > std::numeric_limits<double>::min())
      ++non_zero_terms;
    coeffs.erase(coeffs.begin() + non_zero_terms, coeffs.end());
    
    for (size_t i = 1; i < coeffs.size(); i += 2)
      coeffs[i] *= -1;
    
    if (VERBOSE) {
      cerr << "TERMS USED IN ESTIMATION:" << endl;
      for (size_t i = 0; i < coeffs.size(); ++i)
	cerr << i + 1 << "\t" << coeffs[i] << endl;
      cerr << endl;
    }
    
    const size_t numerator_terms = coeffs.size()/2;
    
    vector<double> numerator, denominator;
    compute_pade_coeffs(coeffs, numerator_terms, numerator, denominator);
    
    if (VERBOSE) {
      cerr << "PADE NUMERATOR COEFFS:" << endl;
      for (size_t i = 0; i < numerator.size(); ++i)
	cerr << i + 1 << "\t" << numerator[i] << endl;
      cerr << endl;
    }
    
    if (VERBOSE) {
      cerr << "PADE DENOMINATOR COEFFS:" << endl;
      for (size_t i = 0; i < denominator.size(); ++i)
	cerr << i + 1 << "\t" << denominator[i] << endl;
      cerr << endl;
    }

    if (VERBOSE)
      cerr << "TOTAL_READS=" << total_reads << endl
	   << "DISTINCT_READS=" << distinct_reads << endl;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    out.setf(std::ios::fixed, std::ios::floatfield);
    using std::setprecision;
    
    for (double t = start_time; t <= end_time; t += time_increment) {
      const double pade_denominator = 
	compute_pade_approx_denominator(t, denominator);
      const double pade_numerator = 
	compute_pade_approx_numerator(t, numerator);
      const double pade_approx = pade_numerator/pade_denominator;
      const double predicted_distinct = distinct_reads + t*pade_approx;
      const double reads_sampled = total_reads*(1.0 + t);
      const double proportion_novel = predicted_distinct/reads_sampled;
      out << std::fixed << setprecision(2) << t << "\t" 
	  << setprecision(0) << reads_sampled << "\t"
	  << setprecision(0) << predicted_distinct << "\t"
	  << std::scientific << setprecision(3) << pade_numerator << "\t"
	  << std::scientific << pade_denominator << "\t"
	  << std::fixed << setprecision(2) 
	  << proportion_novel << endl;
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


