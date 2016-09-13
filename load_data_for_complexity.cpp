/*    Copyright (C) 2014 University of Southern California and
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

#include "load_data_for_complexity.hpp"

#include <queue>
#include <sstream>
#include <unistd.h>


#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "RNG.hpp"

using std::string;
using std::vector;
using std::priority_queue;
using std::min;
using std::endl;
using std::max;
using std::cerr;
using std::tr1::unordered_map;


//////////////////////////////////////////////////////////////////////
// Data imputation
/////////////////////////////////////////////////////////////////////

static bool
update_pe_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                vector<double> &counts_hist,
                                size_t &current_count){
  // check if reads are sorted
  if (curr_gr.same_chrom(prev_gr) &&
      curr_gr.get_start() < prev_gr.get_start() &&
      curr_gr.get_end() < prev_gr.get_end()) {
    return false;
  }
  
  // check if next read is new, and if so update counts_hist to
  // include current_count
  if (!curr_gr.same_chrom(prev_gr) ||
      curr_gr.get_start() != prev_gr.get_start() ||
      curr_gr.get_end() != prev_gr.get_end()) {

    // histogram is too small, resize
    if (counts_hist.size() < current_count + 1)
      counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    current_count = 1;
  }
  else // next read is same, update current_count
    ++current_count;
  
  return true;
}


static void
update_se_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                const string input_file_name,
                                vector<double> &counts_hist,
                                size_t &current_count){
  // check if reads are sorted
  if (curr_gr.same_chrom(prev_gr) &&
      curr_gr.get_start() < prev_gr.get_start())
    throw SMITHLABException("locations unsorted in: " 
                            + input_file_name);

  if (!curr_gr.same_chrom(prev_gr) ||
      curr_gr.get_start() != prev_gr.get_start())
    // next read is new, update counts_hist to include current_count
    {
      // histogram is too small, resize
      if(counts_hist.size() < current_count + 1)
        counts_hist.resize(current_count + 1, 0.0);
      ++counts_hist[current_count];
      current_count = 1;
    }
  else // next read is same, update current_count
    ++current_count;
}



/////comparison function for priority queue/////////////////

/**************** FOR CLARITY BELOW WHEN COMPARING READS *************/
static inline bool
chrom_greater(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_chrom() > b.get_chrom();
}
static inline bool
same_start(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_start() == b.get_start();
}
static inline bool
start_greater(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_start() > b.get_start();
}
static inline bool
end_greater(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_end() > b.get_end();
}
/******************************************************************************/


struct GenomicRegionOrderChecker {
  bool operator()(const GenomicRegion &prev, const GenomicRegion &gr) const {
    return start_check(prev, gr);
  }
  static bool
  start_check(const GenomicRegion &prev, const GenomicRegion &gr) {
    return (chrom_greater(prev, gr)
            || (prev.same_chrom(gr) && start_greater(prev, gr))
            || (prev.same_chrom(gr) && same_start(prev, gr)
                && end_greater(prev, gr)));
  }
};


typedef priority_queue<GenomicRegion, 
                       vector<GenomicRegion>,
                       GenomicRegionOrderChecker> ReadPQ;


static bool
is_ready_to_pop(const ReadPQ &pq, 
                const GenomicRegion &gr,
                const size_t max_width) {
  return !pq.top().same_chrom(gr) || 
    pq.top().get_end() + max_width < gr.get_start();
}


static void
empty_pq(GenomicRegion &curr_gr, GenomicRegion &prev_gr,
         size_t &current_count, vector<double> &counts_hist,
         ReadPQ &read_pq, const string &input_file_name) {
  
  curr_gr = read_pq.top();
  read_pq.pop();
  
  // update counts hist
  const bool UPDATE_SUCCESS =
    update_pe_duplicate_counts_hist(curr_gr, prev_gr, counts_hist,
                                    current_count);
  if (!UPDATE_SUCCESS) {
    std::ostringstream oss;
    oss << "reads unsorted in: " << input_file_name << "\n"
        << "prev = \t" << prev_gr << "\n"
        << "curr = \t" << curr_gr << "\n"
        << "Increase seg_len if in paired end mode";
    throw SMITHLABException(oss.str());
  }
  prev_gr = curr_gr;
}


/*
 * This code is used to deal with read data in BAM format.
 */
#ifdef HAVE_SAMTOOLS
// switching dependency on bamtools to samtools
#include <SAM.hpp>


size_t
load_counts_BAM_se(const string &input_file_name, 
                   vector<double> &counts_hist) {
  const string mapper = "general";
  SAMReader sam_reader(input_file_name, mapper);
  if(!(sam_reader.is_good()))
    throw SMITHLABException("problem opening input file " 
                            + input_file_name);

  SAMRecord samr;
  sam_reader >> samr;
  size_t n_reads = 1;
  // resize vals_hist, make sure it starts out empty
  counts_hist.clear();
  counts_hist.resize(2, 0.0);
  size_t current_count = 1;

  MappedRead prev_mr, curr_mr;
  prev_mr = samr.mr;

  while (sam_reader >> samr, sam_reader.is_good()) {
    // only convert mapped and primary reads
    if (samr.is_primary && samr.is_mapped) {

      // ignore unmapped reads & secondary alignments
      if (!(samr.is_mapping_paired) ||
          (samr.is_mapping_paired && samr.is_Trich)){
        //only count unpaired reads or the left mate of paired reads
        
        curr_mr = samr.mr;
        update_se_duplicate_counts_hist(curr_mr.r, prev_mr.r, 
                                        input_file_name,
                                        counts_hist, 
                                        current_count);
        
        // update number of reads and prev read
        ++n_reads;
        prev_mr = samr.mr;
      }
    }
  }
  
  // to account for the last read compared to the one before it.
  if (counts_hist.size() < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];
  
  return n_reads;
}

/********Below are functions for merging pair-end reads********/

static bool
merge_mates(const size_t suffix_len, const size_t range,
            const GenomicRegion &one, const GenomicRegion &two,
            GenomicRegion &merged, int &len) {

  assert(one.same_chrom(two));
  const size_t read_start = min(one.get_start(), two.get_start());
  const size_t read_end = max(one.get_end(), two.get_end());

  len = read_end - read_start;

  if (len < 0) {
    // cerr << one << endl;
    // cerr << two << endl;
    return false;
  }

  merged = one;
  merged.set_start(read_start);
  merged.set_end(read_end);
  merged.set_score(one.get_score() + two.get_score());

  const string name(one.get_name());
  merged.set_name("FRAG:" + name.substr(0, name.size() - suffix_len));

  return true;
}

inline static bool
same_read(const size_t suffix_len,
          const MappedRead &a, const MappedRead &b) {
  const string sa(a.r.get_name());
  const string sb(b.r.get_name());
  bool SAME_NAME = false;
  if(sa == sb)
    SAME_NAME = true;
  return (SAME_NAME && a.r.same_chrom(b.r));
}

// return true if the genomic region is null
static inline bool
GenomicRegionIsNull(const GenomicRegion &gr){
  GenomicRegion null_gr;
  if(gr == null_gr)
    return true;

  return false;
}


static void
empty_pq(GenomicRegion &prev_gr,
         priority_queue<GenomicRegion,
                        vector<GenomicRegion>,
                        GenomicRegionOrderChecker> &read_pq,
           const string &input_file_name,
         vector<double> &counts_hist,
         size_t &current_count) {

  GenomicRegion curr_gr = read_pq.top();
  read_pq.pop();

  // check if reads are sorted
  if (curr_gr.same_chrom(prev_gr) &&
      curr_gr.get_start() < prev_gr.get_start()
      && curr_gr.get_end() < prev_gr.get_end()) {
    std::ostringstream oss;
    oss << "reads unsorted in: " << input_file_name << "\n"
        << "prev = \t" << prev_gr << "\n"
        << "curr = \t" << curr_gr << "\n"
        << "Increase seg_len if in paired end mode";
    throw SMITHLABException(oss.str());
  }

  if (GenomicRegionIsNull(prev_gr))
    current_count = 1;
  else {
    std::ostringstream oss;
    bool UPDATE_HIST =
      update_pe_duplicate_counts_hist(curr_gr, prev_gr, 
                                      counts_hist, current_count);
    if (!UPDATE_HIST) {
      oss << "locations unsorted in: " << input_file_name << "\n"
          << "prev = \t" << prev_gr << "\n"
          << "curr = \t" << curr_gr << "\n";
      throw SMITHLABException(oss.str());
    }
  }
  prev_gr = curr_gr;
}


size_t
load_counts_BAM_pe(const bool VERBOSE,
                   const string &input_file_name,
                   const size_t MAX_SEGMENT_LENGTH,
                   const size_t MAX_READS_TO_HOLD,
                   size_t &n_paired,
                   size_t &n_mates,
                   vector<double> &counts_hist) {
  
  const string mapper = "general";
  SAMReader sam_reader(input_file_name, mapper);

  // check sam_reader
  if(!(sam_reader.is_good()))
    throw SMITHLABException("problem opening input file " + input_file_name);
  
  SAMRecord samr;
  // resize vals_hist, make sure it starts out empty
  counts_hist.clear();
  counts_hist.resize(2, 0.0);
  size_t current_count = 0;
  size_t suffix_len = 0;
  n_paired = 0;
  n_mates = 0;
  size_t n_unpaired = 0;
  size_t progress_step = 1000000;

  GenomicRegion prev_gr;

  std::priority_queue<GenomicRegion, vector<GenomicRegion>,
                      GenomicRegionOrderChecker> read_pq;
  
  unordered_map<string, SAMRecord> dangling_mates;
  
  while ((sam_reader >> samr, sam_reader.is_good())) {
    
    // only convert mapped and primary reads
    if (samr.is_primary && samr.is_mapped) {
      ++n_mates;
      
      // deal with paired-end stuff
      if (samr.is_mapping_paired) {
        
        const size_t name_len = samr.mr.r.get_name().size() - suffix_len;
        const string read_name(samr.mr.r.get_name().substr(0, name_len));
        
        if (dangling_mates.find(read_name) != dangling_mates.end()) {
          // other end is in dangling mates, merge the two mates
          if(same_read(suffix_len, samr.mr, 
                       dangling_mates[read_name].mr)) {
            if (samr.is_Trich)
              std::swap(samr, dangling_mates[read_name]);
            GenomicRegion merged;
            int len = 0;
            const bool MERGE_SUCCESS =
              merge_mates(suffix_len, MAX_SEGMENT_LENGTH,
                          dangling_mates[read_name].mr.r, samr.mr.r,
                          merged, len);
            // merge success!
            if (MERGE_SUCCESS && len >= 0 &&
                len <= static_cast<int>(MAX_SEGMENT_LENGTH)) {
              read_pq.push(merged);
              ++n_paired;
            }
            else {
              // informative error message!
              if (VERBOSE) {
                cerr << "problem merging read "
                     << read_name << ", splitting read" << endl
                     << samr.mr << "\t" << samr.is_mapping_paired << endl
                     << dangling_mates[read_name].mr << "\t" 
		     << dangling_mates[read_name].is_mapping_paired << endl
                     << "To merge, set max segement "
                     << "length (seg_len) higher." << endl;
              }
              read_pq.push(samr.mr.r);
              read_pq.push(dangling_mates[read_name].mr.r);
              n_unpaired += 2;
            }
            dangling_mates.erase(read_name);
          }
          else {
            read_pq.push(samr.mr.r);
            read_pq.push(dangling_mates[read_name].mr.r);
            dangling_mates.erase(read_name);
            n_unpaired += 2;
          }
        }
        else // didn't find read in dangling_mates, store for later
          dangling_mates[read_name] = samr;
      }
      else {
        read_pq.push(samr.mr.r);
        ++n_unpaired;
      }
      
      
      // dangling mates is too large, flush dangling_mates of reads
      // on different chroms and too far away
      if (dangling_mates.size() > MAX_READS_TO_HOLD) {
        unordered_map<string, SAMRecord> tmp;
        for (unordered_map<string, SAMRecord>::iterator itr =
               dangling_mates.begin(); itr != dangling_mates.end(); ++itr) {
          if (itr->second.mr.r.get_chrom() != samr.mr.r.get_chrom()
              || (itr->second.mr.r.get_chrom() == samr.mr.r.get_chrom()
                  && itr->second.mr.r.get_end() 
                  + MAX_SEGMENT_LENGTH < samr.mr.r.get_start())) {
            if(itr->second.seg_len >= 0) {
              read_pq.push(itr->second.mr.r);
              ++n_unpaired;
            }
          }
          else tmp[itr->first] = itr->second;
        }
        std::swap(tmp, dangling_mates);
        tmp.clear();
      }

      
      // now empty the priority queue
      if (!(read_pq.empty()) &&
          is_ready_to_pop(read_pq, samr.mr.r, MAX_SEGMENT_LENGTH)) {
        //begin emptying priority queue
        while (!(read_pq.empty()) &&
               is_ready_to_pop(read_pq, samr.mr.r, MAX_SEGMENT_LENGTH)) {
          empty_pq(prev_gr, read_pq, input_file_name, counts_hist, current_count);
        }
      }
      
      if (VERBOSE && n_mates % progress_step == 0)
        cerr << "Processed " << n_mates << " records" << endl;
    }
  }

  // empty dangling mates of any excess reads
  while (!dangling_mates.empty()) {
    read_pq.push(dangling_mates.begin()->second.mr.r);
    dangling_mates.erase(dangling_mates.begin());
    ++n_unpaired;
  }
  
  //final iteration
  while(!read_pq.empty())
    empty_pq(prev_gr, read_pq, input_file_name, counts_hist, current_count);
  
  if (counts_hist.size() < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  
  ++counts_hist[current_count];

  assert((read_pq.empty()));
  
  size_t n_reads = n_unpaired + n_paired;
  
  if (VERBOSE)
    cerr << "paired = " << n_paired << endl
         << "unpaired = " << n_unpaired << endl;

  return n_reads;
}

#endif


/* this code is for BED file input */

size_t
load_counts_BED_se(const string input_file_name, 
                   vector<double> &counts_hist) {
  // resize vals_hist
  counts_hist.clear();
  counts_hist.resize(2, 0.0);

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " + input_file_name);
  
  GenomicRegion curr_gr, prev_gr;
  if (!(in >> prev_gr))
    throw SMITHLABException("problem opening file: " + input_file_name);
  
  size_t n_reads = 1;
  size_t current_count = 1;
  while (in >> curr_gr) {
    update_se_duplicate_counts_hist(curr_gr, prev_gr, input_file_name,
                                    counts_hist, current_count);
    ++n_reads;
    prev_gr.swap(curr_gr);
  }
  
  // to account for the last read compared to the one before it.
  if(counts_hist.size() < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];
  
  return n_reads;
}


size_t
load_counts_BED_pe(const string input_file_name, 
                   vector<double> &counts_hist) {

  // resize vals_hist
  counts_hist.clear();
  counts_hist.resize(2, 0.0);

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " 
                            + input_file_name);

  GenomicRegion curr_gr, prev_gr;
  if (!(in >> prev_gr))
    throw SMITHLABException("problem opening file: " 
                            + input_file_name);

  size_t n_reads = 1;
  size_t current_count = 1;

  //read in file and compare each gr with the one before it
  while (in >> curr_gr) {
    const bool UPDATE_SUCCESS =
      update_pe_duplicate_counts_hist(curr_gr, prev_gr,
                                      counts_hist, current_count);
    if (!UPDATE_SUCCESS)
      throw SMITHLABException("reads unsorted in " + input_file_name);
    
    ++n_reads;
    prev_gr.swap(curr_gr);
  }

  if (counts_hist.size() < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  
  // to account for the last read compared to the one before it.
  ++counts_hist[current_count];
  
  return n_reads;

}

/* text file input */
size_t
load_counts(const string &input_file_name, vector<double> &counts_hist) {
  
  std::ifstream in(input_file_name.c_str());
  if (!in) // if file doesn't open
    throw SMITHLABException("problem opening file: " 
                            + input_file_name);
  
  size_t n_reads = 0;
  while(!in.eof()){
    string buffer;
    getline(in, buffer);
    
    std::istringstream iss(buffer);
    if (iss.good()) {
      double val;
      iss >> val;
      if (val > 0) {
        const size_t count = static_cast<size_t>(val);
        // histogram is too small, resize
        if (counts_hist.size() < count + 1)
          counts_hist.resize(count + 1, 0.0);
        ++counts_hist[count];
        n_reads += count;
      }
      else if (val != 0)
        throw SMITHLABException("problem reading file at line " 
                                + toa(n_reads + 1));
    }
    in.peek();
  }
  return n_reads;
}


//returns number of reads from file containing counts histogram
size_t
load_histogram(const string &filename, vector<double> &counts_hist) {
  
  counts_hist.clear();
  
  std::ifstream in(filename.c_str());
  if (!in) //if file doesn't open
    throw SMITHLABException("could not open histogram: " + filename);
  
  size_t n_reads = 0;
  size_t line_count = 0ul, prev_read_count = 0ul;
  string buffer;
  while (getline(in, buffer)) {
    ++line_count;
    size_t read_count = 0ul;
    double frequency = 0.0;
    std::istringstream is(buffer);
    // error reading input
    if (!(is >> read_count >> frequency))
      throw SMITHLABException("bad histogram line format:\n" + buffer + "\n" +
                              "(line " + toa(line_count) + ")");

    // histogram is out of order?
    if (read_count < prev_read_count)
      throw SMITHLABException("bad line order in file " + filename + "\n" +
                              "(line " + toa(line_count) + ")");
    counts_hist.resize(read_count + 1, 0.0);
    counts_hist[read_count] = frequency;
    prev_read_count = read_count;
    n_reads += static_cast<size_t>(read_count*frequency);
  }

  return n_reads;
}


/////////////////////////////////////////////////////////
// Loading coverage counts
////////////////////////////////////////////////////////


// probabilistically split genomic regions into mutiple
// genomic regions of width equal to bin_size
static void
SplitGenomicRegion(const GenomicRegion &inputGR,
                   Runif &runif, const size_t bin_size,
                   vector<GenomicRegion> &outputGRs){
  
  outputGRs.clear();
  GenomicRegion gr(inputGR);

  double frac = static_cast<double>(gr.get_start() % bin_size)/bin_size;
  const size_t width = gr.get_width();
  
  // ADS: this seems like a bunch of duplicated code just for a single
  // function difference
  if (runif.runif(0.0, 1.0) > frac) {
    gr.set_start(std::floor(static_cast<double>(gr.get_start())/
                            bin_size)*bin_size);
    gr.set_end(gr.get_start() + width);
  }
  else {
    gr.set_start(std::ceil(static_cast<double>(gr.get_start())/
                           bin_size)*bin_size);
    gr.set_end(gr.get_start() + width);
  }

  for(size_t i = 0; i < gr.get_width(); i += bin_size){

    const size_t curr_start = gr.get_start() + i;
    const size_t curr_end 
      = std::min(gr.get_end(), curr_start + bin_size);
    frac = static_cast<double>(curr_end - curr_start)/bin_size;

    if(runif.runif(0.0, 1.0) <= frac){
      GenomicRegion binned_gr(gr.get_chrom(), curr_start, 
                              curr_start + bin_size,
                              gr.get_name(), gr.get_score(),
                              gr.get_strand());

      outputGRs.push_back(binned_gr);
    }
  }
}


// split a mapped read into multiple genomic regions
// based on the number of bases in each
static void
SplitMappedRead(const bool VERBOSE,
                const MappedRead &inputMR,
                Runif &runif,
                const size_t bin_size,
                vector<GenomicRegion> &outputGRs){
  
  outputGRs.clear();

  size_t covered_bases = 0;
  size_t read_iterator = inputMR.r.get_start();
  size_t seq_iterator = 0;
  size_t total_covered_bases = 0;
  
  while (seq_iterator < inputMR.seq.size()) {
    if (inputMR.seq[seq_iterator] != 'N')
      covered_bases++;
    
    // if we reach the end of a bin, probabilistically create a binned read
    // with probability proportional to the number of covered bases
    if (read_iterator % bin_size == bin_size - 1) {
      const double frac = static_cast<double>(covered_bases)/bin_size;
      if (runif.runif(0.0, 1.0) <= frac) {
        const size_t curr_start = read_iterator - (read_iterator % bin_size);
        const size_t curr_end = curr_start + bin_size;
        const GenomicRegion binned_gr(inputMR.r.get_chrom(), curr_start, curr_end,
                                      inputMR.r.get_name(), inputMR.r.get_score(),
                                      inputMR.r.get_strand());
        outputGRs.push_back(binned_gr);
      }
      total_covered_bases += covered_bases;
      covered_bases = 0;
    }
    seq_iterator++;
    read_iterator++;
  }

  const double frac = static_cast<double>(covered_bases)/bin_size;
  if (runif.runif(0.0, 1.0) <= frac) {
    const size_t curr_start = read_iterator - (read_iterator % bin_size);
    const size_t curr_end = curr_start + bin_size;
    const GenomicRegion binned_gr(inputMR.r.get_chrom(), curr_start, curr_end,
                                  inputMR.r.get_name(), inputMR.r.get_score(),
                                  inputMR.r.get_strand());
    outputGRs.push_back(binned_gr);
  }
}


size_t
load_coverage_counts_MR(const bool VERBOSE,
                        const string input_file_name,
                        const size_t bin_size,
                        const size_t max_width,
                        vector<double> &coverage_hist,
                        const unsigned int seed) {

  if(VERBOSE){
      cerr << "Setting the random seed to " << seed << endl;
  }
  srand(seed);
  Runif runif(rand());

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " + input_file_name);
  
  MappedRead mr;
  if (!(in >> mr))
    throw SMITHLABException("problem reading from: " + input_file_name);
  
  // initialize prioirty queue to reorder the split reads
  ReadPQ PQ;

  size_t n_reads = 0;
  size_t n_bins = 0;
  GenomicRegion curr_gr, prev_gr;
  size_t current_count = 1;
  
  do {
    
    if (mr.r.get_width() > max_width)
      throw SMITHLABException("Encountered read of width " + 
                              toa(mr.r.get_width()) +
                              "max_width set too small");
    
    vector<GenomicRegion> splitGRs;
    SplitMappedRead(VERBOSE, mr, runif, bin_size, splitGRs);

    n_reads++;
    n_bins += splitGRs.size();

    // add split Genomic Regions to the priority queue
    for (size_t i = 0; i < splitGRs.size(); i++)
      PQ.push(splitGRs[i]);

    // remove Genomic Regions from the priority queue
    if (splitGRs.size() > 0)
      while (!PQ.empty() && is_ready_to_pop(PQ, splitGRs.back(), max_width))
        empty_pq(curr_gr, prev_gr, current_count, 
                 coverage_hist, PQ, input_file_name);
    
  } 
  while (in >> mr);

  // done adding reads, now spit the rest out
  while (!PQ.empty())
    empty_pq(curr_gr, prev_gr, current_count,
             coverage_hist, PQ, input_file_name);
  
  return n_reads;
}


size_t
load_coverage_counts_GR(const string input_file_name,
                        const size_t bin_size,
                        const size_t max_width,
                        vector<double> &coverage_hist,
                        const unsigned int seed) {

  cerr << "Setting the random seed to " << seed << endl;
  srand(seed);
  Runif runif(rand());

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw "problem opening file: " + input_file_name;

  GenomicRegion inputGR;
  if (!(in >> inputGR))
    throw "problem reading from: " + input_file_name;

  // initialize prioirty queue to reorder the split reads
  ReadPQ PQ;

  // prev and current Genomic Regions to compare
  GenomicRegion curr_gr, prev_gr;
  size_t n_reads = 0;
  size_t current_count = 1;

  do {
    
    vector<GenomicRegion> splitGRs;
    SplitGenomicRegion(inputGR, runif, bin_size, splitGRs);
    
    // add split Genomic Regions to the priority queue
    for(size_t i = 0; i < splitGRs.size(); i++)
      PQ.push(splitGRs[i]);
    
    if (splitGRs.size() > 0) {
      // remove Genomic Regions from the priority queue
      while (!PQ.empty() && is_ready_to_pop(PQ, splitGRs.back(), max_width))
        empty_pq(curr_gr, prev_gr, current_count,
                 coverage_hist, PQ, input_file_name);
    }
    n_reads++;
  } 
  while (in >> inputGR);
  
  // done adding reads, now spit the rest out
  while (!PQ.empty())
    empty_pq(curr_gr, prev_gr, current_count,
             coverage_hist, PQ, input_file_name);
  
  return n_reads;
}
