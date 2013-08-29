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
#include <sys/types.h>
#include <unistd.h>


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
		   const BamAlignment &ba) { //passing in paramaters: reference to an unordered map "chrom_lookup", reference to BamAlignment "ba"
	const unordered_map<size_t, string>::const_iterator //constant_iterator "the_chrom" of an unordered map (this is a forward iterator, accesses elements in the direction that goes from beginning to end)
    the_chrom(chrom_lookup.find(ba.RefID)); 	//searches the container for the element that responds to the key "ba.RefID", which is the ID number for the reference sequence
	if (the_chrom == chrom_lookup.end()) //if the resulting element is the same as the "past the end" position in the unordered map chrom_lookup, then 
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID)); //toa converts the RefID integer to string. so if the element is not in the unordered map, prints "no chrom with id REFIDNUMBER" 

  const string chrom = the_chrom->second; // string "chrom" is the element in "the_chrom" that the key points to
  const size_t start = ba.Position; //start of type size_t (unsigned type able to represent the size of any object in bytes) that is set equal to the position (0-based) where the alignment starts
  const size_t end = start + ba.Length; //end is set equal to the position where the alignment starts plus the length of the sequence

  return SimpleGenomicRegion(chrom, start, end); //returns the information for this sequence
}


// same as above, but for paired end reads
static GenomicRegion
BamToGenomicRegion(const unordered_map<size_t, string> &chrom_lookup, 
		   const BamAlignment &ba){ //passing in paramaters: reference to an unordered map "chrom_lookup", reference to BamAlignment "ba"

	const unordered_map<size_t, string>::const_iterator //constant_iterator "the_chrom" of an unordered map (this is a forward iterator, accesses elements in the direction that goes from beginning to end)
    the_chrom(chrom_lookup.find(ba.RefID)); //searches the container for the element that responds to the key "ba.RefID", which is the ID number for the reference sequence
  if (the_chrom == chrom_lookup.end()) //if the resulting element is the same as the "past the end" position in the unordered map chrom_lookup, then 
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID)); //toa converts the RefID integer to string. so if the element is not in the unordered map, prints "no chrom with id REFIDNUMBER"
  
  const string chrom = the_chrom->second; // string "chrom" is the element in "the_chrom" that the key points to
  const size_t start = ba.Position; //start of type size_t (unsigned type able to represent the size of any object in bytes) that is set equal to the position (0-based) where the alignment starts
	const size_t end = ba.Position + ba.InsertSize; //end is set equal to the position where the alignment starts plus the size of the insert

  return GenomicRegion(chrom, start, end); //returns the information for this sequence

}



// loads single end BAM file that returns the number of reads in type size_t (bytes)
static size_t
load_values_BAM_se(const string &input_file_name, vector<double> &values) { //pass in parameters: input file name, vector to hold the counts (unordered)

  BamReader reader; //create a reader
  reader.Open(input_file_name); //opens the BAM file

  // Get header and reference
  string header = reader.GetHeaderText(); 
  RefVector refs = reader.GetReferenceData(); //vector "refs" to hold reference sequence entries

  unordered_map<size_t, string> chrom_lookup; //create unordered map "chrom_lookup"
  for (size_t i = 0; i < refs.size(); ++i) //create variable "i" of type size_t, set to zero. when i < size of "refs" vector, 
    chrom_lookup[i] = refs[i].RefName; //set element in position "i" of chrom_lookup equal to the Name of the reference sequence in position "i" of "refs" vector

  size_t n_reads = 0; //create variable "n_reads" of type size_t, set equal to 0. 
  values.push_back(1.0); //take the vector holding the unordered counts, and add the value of 1.0 to the end of the vector

  SimpleGenomicRegion prev; // create object "prev" of class SimpleGenomicRegion (used for single end)
  BamAlignment bam; // create "bam" object of class BamAlignment, which provides methods to query/modify BAM alignment data fields
  while (reader.GetNextAlignment(bam)) { //while the reader reads through the alignment record from the input file (and outputs it to the alignment destination in "bam") 
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
     //only count unpaired reads or the left mate of paired reads
      if(!(bam.IsPaired()) || 
	 (bam.IsFirstMate())){

	SimpleGenomicRegion r(BamToSimpleGenomicRegion(chrom_lookup, bam)); // create object "r" of class SimpleGenomicRegion, populate it with information about ID, start, and end position of the different sequences
	if (r.same_chrom(prev) && r.get_start() < prev.get_start()) // if sequence in "r" is on the same chromosome as sequence in "prev" AND the start position of sequence in "r" is before the start position of that in "prev", then
	  throw SMITHLABException("locations unsorted in: " + input_file_name); // output message that the BAM file is unsorted. 
    
	if (!r.same_chrom(prev) || r.get_start() != prev.get_start()) // if sequence in "r" is not on the same chromosome as sequence in "prev" OR the start positions are not the same, then 
	  values.push_back(1.0); // add a new value of "1.0" to the end of the vector containing the unordered counts
	else values.back()++; //else (if its in the same chromosome and same start position, we are assuming this read corresponds to the same molecule) ++ to the last element in the vector
	++n_reads; //number of reads increases by 1
	prev.swap(r); // fill object "prev" with sequence information from object "r", and swap information from "prev" into "r" (the first time through the loop, prev is empty, so all the if statements are skipped and this is filled first)
      }
    }
  }
	
	//so in this loop, each sequence is being compared to the sequence before it to determine whether or not it is a read corresponding to the same molecule or not. because the files are sorted, it is sufficient to compare a sequence with the one directly preceding it. 
	//same goes for any code after this that is similar to what is written above
  reader.Close(); //close the reader

  return n_reads; // this function returns the number of reads in the BAM file
}



//loads paired end BAM file and returns the number of reads in type size_t (bytes)
static size_t
load_values_BAM_pe(const string &input_file_name, vector<double> &values) { //pass in parameters: input file name, vector to hold the counts (unordered)

  BamReader reader; //create a reader
  reader.Open(input_file_name); //opens the BAM file

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData(); //vector "refs" to hold reference sequence entries

  unordered_map<size_t, string> chrom_lookup; //create unordered map "chrom_lookup"
  for (size_t i = 0; i < refs.size(); ++i)//create variable "i" of type size_t, set to zero. when i < size of "refs" vector, 
    chrom_lookup[i] = refs[i].RefName;//set element in position "i" of chrom_lookup equal to the Name of the reference sequence in position "i" of "refs" vector

  size_t n_reads = 0; //create variable "n_reads" of type size_t, set equal to 0. 
  values.push_back(1.0); //take the vector holding the unordered counts, and add the value of 1.0 to the end of the vector

  GenomicRegion prev; // create object "prev" of class GenomicRegion (used for paired end)
  BamAlignment bam;  // create "bam" object of class BamAlignment, which provides methods to query/modify BAM alignment data fields
  while (reader.GetNextAlignment(bam)) { //while the reader reads through the alignment record from the input file (and outputs it to the alignment destination in "bam") 
    // ignore unmapped reads & secondary alignments
    if(bam.IsMapped() && bam.IsPrimaryAlignment()){ 
      // ignore reads that do not map concoordantly
      if(bam.IsPaired() && bam.IsProperPair() && bam.IsFirstMate()){ //count only the left mate of concordantly mapped paired end reads
	GenomicRegion r(BamToGenomicRegion(chrom_lookup, bam)); // create object "r" of class GenomicRegion, populate it with information about ID, start, and end position of the different sequences
	if (r.same_chrom(prev) && r.get_start() < prev.get_start() // if sequence in "r" is on the same chromosome as sequence in "prev" AND the start position of sequence in "r" is before the start position of that in "prev" 
	    && r.get_end() < prev.get_end()) //AND if the end position of the sequence in "r" is before the end position of the sequence in "prev", then 
	    throw SMITHLABException("locations unsorted in: " + input_file_name); //output message that the BAM file is unsorted 
    
	if (!r.same_chrom(prev) || r.get_start() != prev.get_start() // if sequence in "r" is not on the same chromosome as sequence in "prev" OR the start positions are not the same
	    || r.get_end() != prev.get_end()) //OR the end positions are not the same, then 
	  values.push_back(1.0);  // add a new value of "1.0" to the end of the vector containing the unordered counts
	else values.back()++; //else ++ to the last element in the vector of unordered counts
	++n_reads; //number of reads increases by 1
	prev.swap(r); // fill object "prev" with sequence information from object "r", and swap information from "prev" into "r" (the first time through the loop, prev is empty, so all the if statements are skipped and this is filled first)
      }
    }
  }
  reader.Close(); //close reader
  return n_reads; //this function returns the number of reads in the BAM file 
}
#endif



//loads single end BED file and returns number of reads in type size_t 
static size_t
load_values_BED_se(const string input_file_name, vector<double> &values) { //pass in parameters: input file name, vector to hold counts (unordered)

 std::ifstream in(input_file_name.c_str()); 
 if (!in) // if file does not open
   throw "problem opening file: " + input_file_name; //error message

 SimpleGenomicRegion r, prev; // create two objects "r" and "prev" of class SimpleGenomicRegion
 if (!(in >> prev)) //read in file to object "prev", but if object "prev" cannot read in file
   throw "problem reading from: " + input_file_name; //error message 

 size_t n_reads = 1; //create variable "n_reads" of type size_t, set equal to 1. 
 values.push_back(1.0); //take the vector holding the unordered counts, and add the value of 1.0 to the end of the vector
 while (in >> r) { //while "r" is reading in the next line of BED file
    if (r.same_chrom(prev) && r.get_start() < prev.get_start())  // if sequence in "r" is on the same chromosome as sequence in "prev" AND the start position of sequence in "r" is before the start position of that in "prev"
      throw SMITHLABException("locations unsorted in: " + input_file_name); //out put error message that the BED file is unsorted 
    
   if (!r.same_chrom(prev) || r.get_start() != prev.get_start())   // if sequence in "r" is not on the same chromosome as sequence in "prev" OR the start positions are not the same
	 values.push_back(1.0); // add a new value of "1.0" to the end of the vector containing the unordered counts
   else values.back()++; //else ++ to the last element in the vector of unordered counts
   ++n_reads; //number of reads increases by 1
   prev.swap(r); // fill object "prev" with sequence information from object "r", and swap information from "prev" into "r" (the first time through the loop, prev is empty, so all the if statements are skipped and this is filled first)
 }
 return n_reads; //returns number of reads in the BED file 
	
	//similar to the above functions for BAM files, this compares previous sequence with the current sequence to create the counts vector
}


//same as above function except for paired end.. 
static size_t
load_values_BED_pe(const string input_file_name, vector<double> &values) {

 std::ifstream in(input_file_name.c_str());
 if (!in)
   throw "problem opening file: " + input_file_name;

 GenomicRegion r, prev;
	in.operator>>(prev);
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



// returns number of reads from file containing observed counts
static size_t
load_values(const string input_file_name, vector<double> &values) { //pass in parameters: input file containing counts, vector to hold counts (unordered) 

  std::ifstream in(input_file_name.c_str());
  if (!in) // if file doesn't open 
    throw SMITHLABException("problem opening file: " + input_file_name); //error message

  vector<double> full_values; //create "full_values" vector, used to hold the values which we will fill the "values" vector with later
  size_t n_reads = 0; 
  static const size_t buffer_size = 10000; // Magic!
  while(!in.eof()){ //as long as it is not end of file, 
    char buffer[buffer_size]; //create array "buffer" of type "char" that has size of 10000
    in.getline(buffer, buffer_size); //read in 10000 chars of the input file and store in "buffer" array
    double val = atof(buffer); //takes string in "buffer" and converts to a double
    if(val > 0.0) //if the value is positive, then 
      full_values.push_back(val); //add this to the "full values" vector

    ++n_reads; //number of reads increases by 1
    in.peek();//check for end of file, if not, next character will be extracted next
  }
  in.close(); //close ifstream
  if(full_values.back() == 0) //if the last element of the "full values" vector is 0, then 
    full_values.pop_back();// delete the zero

  values.swap(full_values); //fill "values" vector with the elements in "full values"
  return n_reads; //return the number of reads
}



//returns number of reads from file containing counts histogram
static void
load_histogram(const string &filename, vector<double> &hist) { //pass in parameters: input file containing counts histogram, vector to hold the histogram

  std::ifstream in(filename.c_str());
  if (!in) //if file doesn't open
    throw SMITHLABException("could not open histogram: " + filename); //error message

  size_t line_count = 0ul, prev_read_count = 0ul; //create variables for prev_read_count and line_count set to 0 unsigned long
  string buffer; //create a string 
  while (getline(in, buffer)) { //while file is being read into "buffer" line by line 
    ++line_count; //increase line count by one as file reads
    size_t read_count = 0ul; //create "read_count", set to 0
    double frequency = 0.0; 
    std::istringstream is(buffer);// use istringstream to copy the string given to it in "buffer" 
    if (!(is >> read_count >> frequency)) //if "is" doesn't input into read_count, and frequency, then  
      throw SMITHLABException("bad histogram line format:\n" +
                              buffer + "\n(line " + toa(line_count) + ")"); //error shows along with line count
    if (read_count < prev_read_count) //if current read count (that has been read in from the line in the file) is less than previous read count, then 
      throw SMITHLABException("bad line order in file " +
                              filename + "\n(line " +
                              toa(line_count) + ")"); //line where mistake is gets printed
    hist.resize(read_count, 0ul); //vector resized to size of histogram (read_count) 
    hist.push_back(frequency); //add frequency to the histogram corresponding to the read_count number
    prev_read_count = read_count; //set prev_read_count to equal read_count (helps keep things in order through the loop)
  }
}

//return how many distinct counts there are 
static double
sample_count_distinct(const gsl_rng *rng,
		      const vector<size_t> &full_umis,
		      const size_t sample_size) { //pass in parameters, a random number generator, a vector for UMIs, and the sample size
  vector<size_t> sample_umis(sample_size); //create vector "sample_umis" of size passed in through parameters 
  gsl_ran_choose(rng, (size_t *)&sample_umis.front(), sample_size,
		 (size_t *)&full_umis.front(), full_umis.size(), 
		 sizeof(size_t)); //fills sample_umis vector with however many elements indicated in sample_size, which are each of size "sizeof(size_t)", taken randomly from the elements in full_umis
  double count = 1.0; 
  for (size_t i = 1; i < sample_umis.size(); i++) //create variable i set equal to 1; so long as this is less than the size of sample_umis, then 
    if(sample_umis[i] != sample_umis[i-1]) //if an element in sample_umis is not equal to the element in the previous position, 
      count++; //counts increases by 1

  return count; //return number of distinct counts 
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
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());
    
    if (VERBOSE)
      cerr << "loading mapped locations" << endl;

    vector<double> values;
    if(HIST_INPUT){
      vector<double> counts_hist;
      load_histogram(input_file_name, counts_hist);

      double total_reads = 0.0;
      for(size_t i = 0; i < counts_hist.size(); i++)
	total_reads += i*counts_hist[i];

      if(VERBOSE){
	cerr << "TOTAL READS     = " << total_reads << endl
	     << "DISTINCT READS  = " << accumulate(counts_hist.begin(), 
						   counts_hist.end(), 0.0) << endl
	     << "MAX COUNT       = " << counts_hist.size() - 1 << endl
	     << "COUNTS OF 1     = " << counts_hist[1] << endl;

	cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
	for (size_t i = 0; i < counts_hist.size(); i++)
	  if (counts_hist[i] > 0)
	    cerr << i << '\t' << counts_hist[i] << endl;
	cerr << endl;
      }

      for(size_t i = 1; i < counts_hist.size(); i++)
	for(size_t j = 0; j < counts_hist[i]; j++)
	  values.push_back(static_cast<double>(i));
    }
    else{
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
