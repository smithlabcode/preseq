/*    c_curve: plot a complexity curve by subsamping sequenced reads
 *    and counting UMIs
 *
 *    Copyright (C) 2013 University of Southern California and
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
#include <sys/types.h>
#include <unistd.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <smithlab_os.hpp>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::tr1::unordered_map;

/*
 * This code is used to deal with read data in BAM format.
 */
#ifdef HAVE_BAMTOOLS
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;



//create BamToSimpleGenomicRegion of class SimpleGenomicRegion (in GenomicRegion.hpp)
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


// same as above, but for paired end reads
static GenomicRegion
BamToGenomicRegion(const unordered_map<size_t, string> &chrom_lookup, 
		   const BamAlignment &ba){

  const unordered_map<size_t, string>::const_iterator
    the_chrom(chrom_lookup.find(ba.RefID)); 
  if (the_chrom == chrom_lookup.end()) 
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  
  const string chrom = the_chrom->second;
  const size_t start = ba.Position; 
  const size_t end = ba.Position + ba.InsertSize; 

  return GenomicRegion(chrom, start, end); 
}



// loads single end BAM file that returns the number of reads 
static size_t
load_values_BAM_se(const string &input_file_name, vector<double> &vals_hist) { 

  // resize vals_hist
  vals_hist.clear();
  vals_hist.resize(2, 0.0);

  BamReader reader; 
  reader.Open(input_file_name); 

  // Get header and reference
  string header = reader.GetHeaderText(); 
  RefVector refs = reader.GetReferenceData(); 

  unordered_map<size_t, string> chrom_lookup; 
  for (size_t i = 0; i < refs.size(); ++i) 
    chrom_lookup[i] = refs[i].RefName; 

  // first read goes in prev, count starts at 1
  BamAlignment bam; 
  reader.GetNextAlignment(bam);
  SimpleGenomicRegion prev(BamToSimpleGenomicRegion(chrom_lookup, bam));
  size_t current_count = 1;
  size_t n_reads = 1; 

  while (reader.GetNextAlignment(bam)) { 
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
     //only count unpaired reads or the left mate of paired reads
      if(!(bam.IsPaired()) || (bam.IsFirstMate())){
	SimpleGenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam)); 
	// check if reads are sorted
	if (r.same_chrom(prev) && r.get_start() < prev.get_start()) 
	  throw SMITHLABException("locations unsorted in: " + input_file_name); 
	// consecutive reads are not duplicates, update histogram
	if (!r.same_chrom(prev) || r.get_start() != prev.get_start()){
	  // histogram is too small, resize
	  if(vals_hist.size() < current_count + 1)
	    vals_hist.resize(current_count + 1, 0.0);
	  ++vals_hist[current_count];
	  current_count = 1;
	}
	else
	  ++current_count;

	++n_reads;
	prev.swap(r); 
      }
    }
  }
	

  reader.Close(); 

  return n_reads; 
}



//loads paired end BAM file and returns the number of reads 
static size_t
load_values_BAM_pe(const string &input_file_name, vector<double> &vals_hist) { 

  // resize vals_hist
  vals_hist.clear();
  vals_hist.resize(2, 0.0);

  BamReader reader;
  reader.Open(input_file_name); 

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData(); 

  unordered_map<size_t, string> chrom_lookup; 
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;

  // first read goes in prev, count starts at 1
  BamAlignment bam; 
  reader.GetNextAlignment(bam);
  GenomicRegion prev(BamToGenomicRegion(chrom_lookup, bam));
  size_t current_count = 1;
  size_t n_reads = 1; 

  while (reader.GetNextAlignment(bam)) { 
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
      // ignore reads that do not map concoordantly
      if(bam.IsPaired() && bam.IsProperPair() && bam.IsFirstMate()){ 
	GenomicRegion r(BamToGenomicRegion(chrom_lookup, bam));
	// check if reads are sorted
	if (r.same_chrom(prev) && r.get_start() < prev.get_start() && r.get_end() < prev.get_end()) 
	    throw SMITHLABException("locations unsorted in: " + input_file_name); 
    
	if (!r.same_chrom(prev) || r.get_start() != prev.get_start() || r.get_end() != prev.get_end()) {
	  // histogram is too small, resize
	  if(vals_hist.size() < current_count + 1)
	    vals_hist.resize(current_count + 1, 0.0);
	  ++vals_hist[current_count];
	  current_count = 1;
	}
	else
	  ++current_count;

	++n_reads; 
	prev.swap(r); 
      }
    }
  }
  reader.Close(); 
  return n_reads; 
}
#endif



//loads single end BED file and returns number of reads 
static size_t
load_values_BED_se(const string input_file_name, vector<double> &vals_hist) { 

  // resize vals_hist
  vals_hist.clear();
  vals_hist.resize(2, 0.0);

 std::ifstream in(input_file_name.c_str()); 
 if (!in) // if file does not open
   throw "problem opening file: " + input_file_name; 

 SimpleGenomicRegion r, prev; 
 if (!(in >> prev)) // problem reading
   throw "problem reading from: " + input_file_name; 

 size_t n_reads = 1; 
 size_t current_count = 1;

 while (in >> r) {
   // check if reads are sorted
   if (r.same_chrom(prev) && r.get_start() < prev.get_start()) 
     throw SMITHLABException("locations unsorted in: " + input_file_name);
    
   if (!r.same_chrom(prev) || r.get_start() != prev.get_start()) {  
     // histogram is too small, resize
     if(vals_hist.size() < current_count + 1)
       vals_hist.resize(current_count + 1, 0.0);
     ++vals_hist[current_count];
     current_count = 1;
   }
   else
     ++current_count;

   ++n_reads; 
   prev.swap(r); 
 }

 return n_reads; 
}


//same as above function except for paired end.. 
static size_t
load_values_BED_pe(const string input_file_name, vector<double> &vals_hist) {

  // resize vals_hist
  vals_hist.clear();
  vals_hist.resize(2, 0.0);

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 GenomicRegion r, prev;
 if (!(in >> prev))
   throw "problem reading from: " + input_file_name;

 size_t n_reads = 1;
 size_t current_count = 1;

 while (in >> r) {
   // check if reads are sorted
   if (r.same_chrom(prev) && r.get_start() < prev.get_start() && r.get_end() < prev.get_end())
     throw SMITHLABException("locations unsorted in: " + input_file_name);
    
   if (!r.same_chrom(prev) || r.get_start() != prev.get_start() || r.get_end() != prev.get_end()) {
     // histogram is too small, resize
     if(vals_hist.size() < current_count + 1)
       vals_hist.resize(current_count + 1, 0.0);
     ++vals_hist[current_count];
     current_count = 1;
   }
   else
     ++current_count;

   ++n_reads;
   prev.swap(r);
 }
 return n_reads;
}



// returns number of reads from file containing observed counts
static size_t
load_values(const string input_file_name, vector<double> &vals_hist) { 

  std::ifstream in(input_file_name.c_str());
  if (!in) // if file doesn't open 
    throw SMITHLABException("problem opening file: " + input_file_name); //error message

  vector<double> values; 
  size_t n_reads = 0; 
  static const size_t buffer_size = 10000; // Magic!
  while(!in.eof()){ 
    char buffer[buffer_size]; 
    in.getline(buffer, buffer_size);
    double val = atof(buffer); 
    if(val > 0.0) 
      values.push_back(val); 

    ++n_reads; 
    in.peek();
  }
  in.close(); 

  const size_t max_observed_count = 
    static_cast<size_t>(*std::max_element(values.begin(), values.end())); 
  vector<double> counts_hist(max_observed_count + 1, 0.0); 
  for (size_t i = 0; i < values.size(); ++i)
    ++counts_hist[static_cast<size_t>(values[i])];

  vals_hist.swap(counts_hist); 
  return n_reads; 
}



//returns number of reads from file containing counts histogram
static void
load_histogram(const string &filename, vector<double> &hist) { 
  
  hist.clear();

  std::ifstream in(filename.c_str());
  if (!in) //if file doesn't open
    throw SMITHLABException("could not open histogram: " + filename); 

  size_t line_count = 0ul, prev_read_count = 0ul; 
  string buffer;
  while (getline(in, buffer)) { 
    ++line_count; 
    size_t read_count = 0ul; 
    double frequency = 0.0; 
    std::istringstream is(buffer);
    // error reading input
    if (!(is >> read_count >> frequency)) 
      throw SMITHLABException("bad histogram line format:\n" +
                              buffer + "\n(line " + toa(line_count) + ")"); 
    // histogram is out of order
    if (read_count < prev_read_count) 
      throw SMITHLABException("bad line order in file " +
                              filename + "\n(line " +
                              toa(line_count) + ")"); 
    hist.resize(read_count + 1, 0.0);
    hist[read_count] = frequency; 
    prev_read_count = read_count; 
  }
}

// interpolate complexity
//return how many distinct counts there are in a sample from a full set of UMIs 
static double
sample_count_distinct(const gsl_rng *rng,
		      const vector<size_t> &full_umis,
		      const size_t sample_size) {
  vector<size_t> sample_umis(sample_size); 
  // sample from UMIs w/out replacement
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
    //    size_t lower_limit = 1000000;
    size_t upper_limit = 0;
    size_t step_size = 1000000;

    bool VERBOSE = false;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;


#ifdef HAVE_BAMTOOLS
    bool BAM_FORMAT_INPUT = false;
#endif


    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("c_curve", "plot a complexity curve by subsamping "
			   "sequenced reads and counting UMIs",
			   "<bed-file|bam-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    //        opt_parse.add_opt("lower", 'l', "lower limit for samples", 
    //		      false , lower_limit);
    //	      opt_parse.add_opt("upper", 'u', "upper limit for samples", 
    //		      false , upper_limit);
        opt_parse.add_opt("step", 's', "step size for samples", 
    		      false , step_size);
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
    opt_parse.add_opt("hist", 'H', 
		      "input is a text file containing the observed histogram",
		      false, HIST_INPUT);

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
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); // use default type 
    srand(time(0) + getpid()); //give the random fxn a new seed
    gsl_rng_set(rng, rand()); //initialize random number generator with the seed 
    
    if (VERBOSE)
      cerr << "loading mapped locations" << endl;

    vector<double> counts_hist;
    if(HIST_INPUT) 
      load_histogram(input_file_name, counts_hist); 
    else if(VALS_INPUT)
	load_values(input_file_name, counts_hist); 
#ifdef HAVE_BAMTOOLS
	//if user decides to input BAM files 
    else if (BAM_FORMAT_INPUT && PAIRED_END) 
      load_values_BAM_pe(input_file_name, counts_hist);
    else if(BAM_FORMAT_INPUT) 
      load_values_BAM_se(input_file_name, counts_hist);
#endif
    // paired end bed or mr format
    else if(PAIRED_END) 
	load_values_BED_pe(input_file_name, counts_hist);  
    // default is single end bed or mr format   
    else
      load_values_BED_se(input_file_name, counts_hist); 

    size_t total_reads = 0ul;
    for(size_t i = 0; i < counts_hist.size(); i++)
      total_reads += static_cast<size_t>(counts_hist[i])*i;

    if(VERBOSE){
      cerr << "TOTAL READS     = " << total_reads << endl
	   << "DISTINCT READS  = " << accumulate(counts_hist.begin(), counts_hist.end(), 0.0) << endl
	   << "MAX COUNT       = " << counts_hist.size() - 1 << endl
	   << "COUNTS OF 1     = " << counts_hist[1] << endl;

      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    }

    vector<size_t> full_umis;
    size_t umi = 1;
    for(size_t i = 1; i < counts_hist.size(); i++){
      for(size_t j = 0; j < counts_hist[i]; j++){
	for(size_t k = 0; k < i; k++)
	  full_umis.push_back(umi);
	umi++;
      }
    }
    // sanity check
    assert(full_umis.size() == total_reads);
    
    if (upper_limit == 0)
      upper_limit = full_umis.size(); //set upper limit to equal the number of molecules
    
    // print curve
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
	 
    out << "TOTAL_READS" << "\t" << "DISTINCT_READS" << endl;
    out << 0 << '\t' << 0 << endl;
    for (size_t i = step_size; i <= upper_limit; i += step_size) { 
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
