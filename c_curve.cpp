/*    c_curve: plot a complexity curve by subsamping sequenced reads
 *    and counting UMIs
 *
 *    Copyright (C) 2012 University of Southern California and
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

#include <iomanip>
#include <numeric>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::tr1::unordered_map;


#ifdef HAVE_BAMTOOLS
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;

static SimpleGenomicRegion
BamToSimpleGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
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

static GenomicRegion
BamToGenomicRegion(const unordered_map<size_t, string> &chrom_lookup,
		   const BamAlignment &ba){

  const unordered_map<size_t, string>::const_iterator
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = ba.InsertSize;

  return GenomicRegion(chrom, start, end);

}


static size_t
load_values_BAM_se(const string &input_file_name, vector<double> &values) {

  BamReader reader;
  reader.Open(input_file_name);

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;

  size_t n_reads = 1;
  values.push_back(1.0);

  SimpleGenomicRegion prev;
  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
      if(bam.IsPaired()) // throw error if paired end read found
	throw SMITHLABException("paired end input unexpectedly found in " 
				+ input_file_name);
      else{
	SimpleGenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam));
	if (r.same_chrom(prev) && r.get_start() < prev.get_start())
	  throw SMITHLABException("locations unsorted in: " + input_file_name);
    
	if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
	  values.push_back(1.0);
	else values.back()++;
	++n_reads;
	prev.swap(r);
      }
    }
  }
  reader.Close();
  return n_reads;
}

static size_t
load_values_BAM_pe(const string &input_file_name, vector<double> &values) {

  BamReader reader;
  reader.Open(input_file_name);

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;

  size_t n_reads = 1;
  values.push_back(1.0);

  GenomicRegion prev;
  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
      // ignore reads that do not map concoordantly
      if(bam.IsPaired() && bam.IsProperPair() && bam.IsFirstMate()){
	GenomicRegion r(BamToGenomicRegion(chrom_lookup, bam));
	if (r.same_chrom(prev) && r.get_start() < prev.get_start()
	    && r.get_end() < prev.get_end())
	    throw SMITHLABException("locations unsorted in: " + input_file_name);
    
	if (!r.same_chrom(prev) || r.get_start() != prev.get_start()
	    || r.get_end() != prev.get_end())
	  values.push_back(1.0);
	else values.back()++;
	++n_reads;
	prev.swap(r);
      }
    }
  }
  reader.Close();
  return n_reads;
}
#endif


static size_t
load_values_BED_se(const string input_file_name, vector<double> &values) {

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 SimpleGenomicRegion r, prev;
 if (!(in >> prev))
   throw "problem reading from: " + input_file_name;

 size_t n_reads = 1;
 values.push_back(1.0);
 while (in >> r) {
    if (r.same_chrom(prev) && r.get_start() < prev.get_start())
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    
   if (!r.same_chrom(prev) || r.get_start() != prev.get_start())
     values.push_back(1.0);
   else values.back()++;
   ++n_reads;
   prev.swap(r);
 }
 return n_reads;
}

static size_t
load_values_BED_pe(const string input_file_name, vector<double> &values) {

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 GenomicRegion r, prev;
 if (!(in >> prev))
   throw "problem reading from: " + input_file_name;

 size_t n_reads = 1;
 values.push_back(1.0);
 while (in >> r) {
    if (r.same_chrom(prev) && r.get_start() < prev.get_start() 
	&& r.get_end() < prev.get_end())
      throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!r.same_chrom(prev) || r.get_start() != prev.get_start() 
	|| r.get_end() != prev.get_end())
     values.push_back(1.0);
   else values.back()++;
   ++n_reads;
   prev.swap(r);
 }
 return n_reads;
}


static size_t
load_values(const string input_file_name, vector<double> &values) {

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " + input_file_name);

  vector<double> full_values;
  size_t n_reads = 0;
  static const size_t buffer_size = 10000; // Magic!
  while(!in.eof()){
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    full_values.push_back(atof(buffer));
    if(full_values.back() <= 0.0){
      cerr << "INVALID INPUT\t" << buffer << endl;
      throw SMITHLABException("ERROR IN INPUT");
    }
    ++n_reads;
    in.peek();
  }
  in.close();
  if(full_values.back() == 0)
    full_values.pop_back();

  values.swap(full_values);
  return n_reads;
}


static double
sample_count_distinct(const gsl_rng *rng,
		      const vector<size_t> &full_umis,
		      const size_t sample_size) {
  vector<size_t> sample_umis(sample_size);
  gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
		 (size_t *)&full_umis.front(), full_umis.size(), 
		 sizeof(size_t));
  double count = 1.0;
  for (size_t i = 1; i < sample_umis.size(); i++)
    if(sample_umis[i] != sample_umis[i-1])
      count++;

  return count;
}


int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    size_t lower_limit = 1000000;
    size_t upper_limit = 0;
    size_t step_size = 1000000;

    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;

#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif


    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("c_curve", "plot a complexity curve by subsamping "
			   "sequenced reads and counting UMIs",
			   "<bed-file|bam-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    //    opt_parse.add_opt("lower", 'l', "lower limit for samples", 
    //		      false , lower_limit);
    //    opt_parse.add_opt("upper", 'u', "upper limit for samples", 
    //		      false , upper_limit);
    //    opt_parse.add_opt("step", 's', "step size for samples", 
    //		      false , step_size);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false , VERBOSE);
#ifdef HAVE_BAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format", 
		      false , BAM_FORMAT_INPUT);
#endif
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
		      false, PAIRED_END);
    opt_parse.add_opt("vals", 'V', 
		      "input is a text file containing only the observed counts",
		      false, VALS_INPUT);

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
    
    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());
    
    if (VERBOSE)
      cerr << "loading mapped locations" << endl;

    vector<double> values;
    if(VALS_INPUT)
      load_values(input_file_name, values);
#ifdef HAVE_BAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END)
      load_values_BAM_pe(input_file_name, values);
    else if(BAM_FORMAT_INPUT)
      load_values_BAM_se(input_file_name, values);
#endif
    else if(PAIRED_END)
      load_values_BED_pe(input_file_name, values);  
    else
      load_values_BED_se(input_file_name, values);

    if(VERBOSE){
      const size_t max_observed_count = 
	static_cast<size_t>(*std::max_element(values.begin(), values.end()));
      vector<double> counts_hist(max_observed_count + 1, 0.0);
      for (size_t i = 0; i < values.size(); ++i)
	++counts_hist[static_cast<size_t>(values[i])];
    
      cerr << "TOTAL READS     = " << accumulate(values.begin(), values.end(), 0.0) << endl
	   << "DISTINCT READS  = " << values.size() << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_hist[1] << endl;

      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    }

    vector<size_t> full_umis;
    for (size_t i = 0; i < values.size(); i++)
      for (size_t j = 0; j < values[i]; j++)
	full_umis.push_back(i+1);
    
    if (upper_limit == 0)
      upper_limit = full_umis.size();
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    out << "total_reads" << "\t" << "distinct_reads" << endl;
    for (size_t i = lower_limit; i <= upper_limit; i += step_size) {
      if (VERBOSE)
	cerr << "sample size: " << i << endl;
      out << i << "\t" << sample_count_distinct(rng, full_umis, i) << endl;
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
