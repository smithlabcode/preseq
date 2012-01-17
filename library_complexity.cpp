/*    library_complexity:
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Timothy Daley
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


#include "pade_approximant.hpp"
#include "continued_fraction.hpp"
#include "library_size_estimates.hpp"

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>

#include <fstream>
#include <numeric>
#include <vector>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::ostream;

using std::setw;
using std::fixed;
using std::setprecision;
using std::tr1::unordered_map;

using smithlab::log_sum_log_vec;

#ifdef HAVE_BAMTOOLS
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;

static SimpleGenomicRegion
BamAlignmentToSimpleGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
				  const BamAlignment &ba) {
  const unordered_map<size_t, string>::const_iterator 
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  
  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = start + ba.Length;
  return SimpleGenomicRegion(chrom, start, end);
}


static void
ReadBAMFormatInput(const string &infile, vector<SimpleGenomicRegion> &read_locations) {
  
  BamReader reader;
  reader.Open(infile);
  
  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  
  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;
  
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
    read_locations.push_back(BamAlignmentToSimpleGenomicRegion(chrom_lookup, bam));
  reader.Close();
}
#endif


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
get_counts(const vector<SimpleGenomicRegion> &reads, vector<size_t> &counts) {
  counts.push_back(1);
  for (size_t i = 1; i < reads.size(); ++i)
    if (reads[i] == reads[i - 1]) counts.back()++;
    else counts.push_back(1);
}



int
main(const int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    string stats_outfile;

    size_t max_terms = 100;
    double max_extrapolation = 1e10;
    double step_size = 1e6;
    size_t smoothing_bandwidth = 4;
    double smoothing_decay_factor = 15.0;

    int diagonal = 0;
    
    /* FLAGS */
    bool VERBOSE = false;
    bool SMOOTH_HISTOGRAM = true; //false; 
    // bool LIBRARY_SIZE = false;
    
#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(argv[0], "", "<sorted-bed-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("stats", 'S', "stats output file", 
		      false , stats_outfile);
    opt_parse.add_opt("extrapolation_length",'e',"maximum extrapolation length", 
                      false, max_extrapolation);
    opt_parse.add_opt("step",'s',"step size between extrapolations", 
                      false, step_size);
    opt_parse.add_opt("terms",'t',"maximum number of terms", false, max_terms);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
#ifdef HAVE_BAMTOOLS
    opt_parse.add_opt("bam", 'b', "input is in BAM format", 
		      false , BAM_FORMAT_INPUT);
#endif
    //     opt_parse.add_opt("tol", '\0', "general numerical tolerance",
    // 		      false, tolerance);
    //     opt_parse.add_opt("delta", '\0', "derivative step size",
    //                       false, deriv_delta);
    //     opt_parse.add_opt("smooth",'\0',"smooth histogram (default: no smoothing)",
    //                       false, SMOOTH_HISTOGRAM);
    //     opt_parse.add_opt("bandwidth", '\0', "smoothing bandwidth",
    // 		      false, smoothing_bandwidth);
    //     opt_parse.add_opt("decay", '\0', "smoothing decay factor",
    // 		      false, smoothing_decay_factor);
    //     opt_parse.add_opt("diag", 'd', "diagonal to use",
    // 		      false, diagonal);
    //     opt_parse.add_opt("library_size", '\0', "estimate library size "
    // 		      "(default: estimate distinct)", false, LIBRARY_SIZE);
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
    
    // READ IN THE DATA
    vector<SimpleGenomicRegion> read_locations;
#ifdef HAVE_BAMTOOLS
    if (BAM_FORMAT_INPUT)
      ReadBAMFormatInput(input_file_name, read_locations);
    else 
#endif
      ReadBEDFile(input_file_name, read_locations);
    if (!check_sorted(read_locations))
      throw SMITHLABException("read_locations not sorted");

    // OBTAIN THE COUNTS FOR DISTINCT READS
    vector<size_t> values;
    get_counts(read_locations, values);
    
    // JUST A SANITY CHECK
    const size_t n_reads = read_locations.size();
    const size_t values_sum = accumulate(values.begin(), values.end(), 0ul);
    assert(values_sum == n_reads);
    
    const double max_val = max_extrapolation/static_cast<double>(values_sum);
    const double val_step = step_size/static_cast<double>(values_sum);
    
    const size_t max_observed_count = 
      *std::max_element(values.begin(), values.end());

    // BUILD THE HISTOGRAM
    vector<double> counts_histogram(max_observed_count + 1, 0.0);
    for (size_t i = 0; i < values.size(); ++i)
      ++counts_histogram[values[i]];
    
    const size_t distinct_counts = 
      std::count_if(counts_histogram.begin(), counts_histogram.end(),
		    bind2nd(std::greater<size_t>(), 0));
    
    // initialize ContinuedFractionApproximation 
    if (SMOOTH_HISTOGRAM) 
      smooth_histogram(smoothing_bandwidth, 
		       smoothing_decay_factor, counts_histogram);
    
    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < counts_histogram.size() && 
	   counts_histogram[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
    max_terms = std::min(max_terms, counts_before_first_zero);     
    
    if (VERBOSE)
      cerr << "TOTAL READS     = " << read_locations.size() << endl
	   << "DISTINCT COUNTS = " << distinct_counts << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_histogram[1] << endl
	   << "MAX TERMS       = " << max_terms << endl;
    
    /* Now do the work associated with the continued fraction or the
       Pade table.
     */
    const ContinuedFractionApproximation cfa(diagonal, max_terms, step_size, max_extrapolation);
    const ContinuedFraction cf(cfa.optimal_continued_fraction(counts_histogram));
    
    if (VERBOSE)
      cerr << cf << endl;
    
    /* Here is where we actually compute the estimates and either
       write them to a file or stdout, as specified by the user.
    */
    vector<double> estimates;
    cf.extrapolate_distinct(counts_histogram, max_val, val_step, estimates);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    double val = 0.0;
    for (size_t i = 0; i < estimates.size(); ++i, val += val_step)
      out << std::fixed << std::setprecision(1) 
	  << (val + 1.0)*values_sum << '\t' << estimates[i] << endl;
    
    /* Finally output the bounds on the library size if either verbose
       output is requested or this info was specifically requested.
     */
    if (VERBOSE || !stats_outfile.empty()) {
      std::ofstream stats_of;
      if (!stats_outfile.empty()) stats_of.open(stats_outfile.c_str());
      ostream stats_out(stats_outfile.empty() ? cerr.rdbuf() : stats_of.rdbuf());
      
      const double upper_bound = distinct_counts +
	upperbound_librarysize(counts_histogram, max_terms);
      
      stats_out << "CHAO87_LOWER=" 
		<< chao87_lowerbound_librarysize(counts_histogram) << endl;
      stats_out	<< "CL92_LOWER=" 
		<< cl92_lowerbound_librarysize(counts_histogram) << endl;
      stats_out << "CF_LOWER="
		<< cfa.lowerbound_librarysize(counts_histogram, upper_bound) << endl;
      stats_out << "CF_UPPER=" << upper_bound << endl;
    }
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
