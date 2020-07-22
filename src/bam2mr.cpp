/*    to-mr: a program for converting SAM and BAM format to MappedRead
 *    format.
 *    Currently supported mappers: bsmap, bismark.
 *
 *    Copyright (C) 2009-2012 University of Southern California and
 *                            Andrew D. Smith
 *
 *    Authors: Meng Zhou, Qiang Song, Andrew Smith
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

#include <cmath>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "SAM.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::max;
using std::min;
using std::priority_queue;
using std::tr1::unordered_map;





/********Below are functions for merging pair-end reads********/
static void
fill_overlap(const bool pos_str, const MappedRead &mr, const size_t start, 
             const size_t end, const size_t offset, string &seq, string &scr) {
  const size_t a = pos_str ? (start - mr.r.get_start()) : (mr.r.get_end() - end);
  const size_t b = pos_str ? (end -  mr.r.get_start()) : (mr.r.get_end() - start);
  copy(mr.seq.begin() + a, mr.seq.begin() + b, seq.begin() + offset);
  copy(mr.scr.begin() + a, mr.scr.begin() + b, scr.begin() + offset);
}

static bool
merge_mates(const size_t suffix_len, const size_t range,
            const MappedRead &one, const MappedRead &two,
            MappedRead &merged, int &len) {
  
  const bool pos_str = one.r.pos_strand();
  const size_t overlap_start = max(one.r.get_start(), two.r.get_start());
  const size_t overlap_end = min(one.r.get_end(), two.r.get_end());

  const size_t one_left = pos_str ? 
    one.r.get_start() : max(overlap_end, one.r.get_start());
  const size_t one_right = 
    pos_str ? min(overlap_start, one.r.get_end()) : one.r.get_end();
  
  const size_t two_left = pos_str ? 
    max(overlap_end, two.r.get_start()) : two.r.get_start();
  const size_t two_right = pos_str ? 
    two.r.get_end() : min(overlap_start, two.r.get_end());

  len = pos_str ? (two_right - one_left) : (one_right - two_left);
  
  if(len < 0){
    cerr << one << endl;
    cerr << two << endl;
    return false;
  }
  // assert(len > 0);
  // assert(one_left <= one_right && two_left <= two_right);
  // assert(overlap_start >= overlap_end || static_cast<size_t>(len) == 
  //    ((one_right - one_left) + (two_right - two_left) + (overlap_end - overlap_start)));
  
  string seq(len, 'N');
  string scr(len, 'B');
  if (len >= 0 && len <= static_cast<int>(range)) {
    // lim_one: offset in merged sequence where overlap starts
    const size_t lim_one = one_right - one_left;
    copy(one.seq.begin(), one.seq.begin() + lim_one, seq.begin());
    copy(one.scr.begin(), one.scr.begin() + lim_one, scr.begin());
    
    const size_t lim_two = two_right - two_left;
    copy(two.seq.end() - lim_two, two.seq.end(), seq.end() - lim_two);
    copy(two.scr.end() - lim_two, two.scr.end(), scr.end() - lim_two);
    
    // deal with overlapping part
    if (overlap_start < overlap_end) {
      const size_t one_bads = count(one.seq.begin(), one.seq.end(), 'N');
      const int info_one = one.seq.length() - (one_bads + one.r.get_score());
      
      const size_t two_bads = count(two.seq.begin(), two.seq.end(), 'N');
      const int info_two = two.seq.length() - (two_bads + two.r.get_score());
      
      // use the mate with the most info to fill in the overlap
      if (info_one >= info_two)
	fill_overlap(pos_str, one, overlap_start, overlap_end, lim_one, seq, scr);
      else
	fill_overlap(pos_str, two, overlap_start, overlap_end, lim_one, seq, scr);
    }
  }
  
  merged = one;
  merged.r.set_start(pos_str ? one.r.get_start() : two.r.get_start());
  merged.r.set_end(merged.r.get_start() + len);
  merged.r.set_score(one.r.get_score() + two.r.get_score());
  merged.seq = seq;
  merged.scr = scr;  
  const string name(one.r.get_name());
  merged.r.set_name("FRAG:" + name.substr(0, name.size() - suffix_len));

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



static void
revcomp(MappedRead &mr) {
  // set the strand to the opposite of the current value
  mr.r.set_strand(mr.r.pos_strand() ? '-' : '+');
  // reverse complement the sequence, and reverse the quality scores
  revcomp_inplace(mr.seq);
  std::reverse(mr.scr.begin(), mr.scr.end());
}
/********Above are functions for merging pair-end reads********/


/////comparison function for priority queue/////////////////

/**************** FOR CLARITY BELOW WHEN COMPARING READS *************/
static inline bool
chrom_greater(const MappedRead &a, const MappedRead &b) {
    return a.r.get_chrom() > b.r.get_chrom();
}
static inline bool
same_start(const MappedRead &a, const MappedRead &b) {
    return a.r.get_start() == b.r.get_start();
}
static inline bool
start_greater(const MappedRead &a, const MappedRead &b) {
    return a.r.get_start() > b.r.get_start();
}
static inline bool
end_greater(const MappedRead &a, const MappedRead &b) {
    return a.r.get_end() > b.r.get_end();
}
/******************************************************************************/


struct MappedReadOrderChecker {
    bool operator()(const MappedRead &prev, const MappedRead &gr) const {
        return start_check(prev, gr);
    }
    static bool
    is_ready(const priority_queue<MappedRead, vector<MappedRead>, MappedReadOrderChecker> &pq,
             const MappedRead &mr, const size_t max_width) {
        return !pq.top().r.same_chrom(mr.r) || pq.top().r.get_end() + max_width < mr.r.get_start();
    }
    static bool
    start_check(const MappedRead &prev, const MappedRead &mr) {
        return (chrom_greater(prev, mr)
                || (prev.r.same_chrom(mr.r) && start_greater(prev, mr))
                || (prev.r.same_chrom(mr.r) && same_start(prev, mr) && end_greater(prev, mr)));
    }
};



static void empty_pq(MappedRead &prev_mr,
                     priority_queue<MappedRead, vector<MappedRead>,
                                    MappedReadOrderChecker> &read_pq,
                     const string &input_file_name,
		     std::ostream& out){
    
  MappedRead curr_mr = read_pq.top();
    //	       cerr << "outputting from queue : " << read_pq.top() << endl;
  read_pq.pop();

  /*
    // check if reads are sorted
  if (curr_mr.r.same_chrom(prev_mr.r) &&
      curr_mr.r.get_start() < prev_mr.r.get_start()
      && curr_mr.r.get_end() < prev_mr.r.get_end()){
    cerr << "prev = \t" << prev_mr << endl;
    cerr << "curr = \t" << curr_mr << endl;
    cerr << "Increase seg_len if in paired end mode" << endl;
    throw SMITHLABException("reads unsorted in " + input_file_name);
  }
  */

  out << curr_mr << endl;

  prev_mr = curr_mr;
}


int 
main(int argc, const char **argv) {
  try {
    string outfile;
    string mapper = "general";
    size_t MAX_SEGMENT_LENGTH = 10000;
    size_t suffix_len = 0;
    bool VERBOSE = false;
    size_t MAX_READS_TO_HOLD = 1000000;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
			   "", "<bam file>");
    opt_parse.add_opt("output", 'o', "Name of output file", 
                      false, outfile);
    opt_parse.add_opt("suff", 's', "read name suffix length (default: 0)",
                      false, suffix_len); 
    opt_parse.add_opt("seg_len", 'L', "maximum allowed insert size", 
                      false, MAX_SEGMENT_LENGTH); 
    opt_parse.add_opt("max_reads", 'R', "maximum number of reads to hold for merging",
		      false, MAX_READS_TO_HOLD);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc < 2 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
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
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if (VERBOSE)
    {
      cerr << "Input file: " << mapped_reads_file << endl
           << "Output file: " << (outfile.empty() ? "stdout" : outfile) << endl;
    }

    SAMReader sam_reader(mapped_reads_file, mapper);
    std::tr1::unordered_map<string, SAMRecord> dangling_mates;
   
    const size_t progress_step = 1000000;
    SAMRecord samr;
    MappedRead prev_mr;
    size_t n_mates = 0;

    std::priority_queue<MappedRead, vector<MappedRead>, MappedReadOrderChecker> read_pq;

    while ((sam_reader >> samr, sam_reader.is_good()))
    {
      if(samr.is_primary && samr.is_mapped){
	// only convert mapped and primary reads
	++n_mates;
	if (samr.is_mapping_paired){
	  const string read_name
	    = samr.mr.r.get_name().substr(0, samr.mr.r.get_name().size() - suffix_len);

	  if (dangling_mates.find(read_name) != dangling_mates.end()){
	    // other end is in dangling mates, merge the two mates
	    if(same_read(suffix_len, samr.mr, dangling_mates[read_name].mr)){
	      if (samr.is_Trich) std::swap(samr, dangling_mates[read_name]);

	      revcomp(samr.mr);

	      MappedRead merged;
	      int len = 0;
	      bool MERGE_SUCCESS =
		merge_mates(suffix_len, MAX_SEGMENT_LENGTH,
			    dangling_mates[read_name].mr, samr.mr, merged, len);

	      if (MERGE_SUCCESS && 
		  len >= 0 && 
		  len <= static_cast<int>(MAX_SEGMENT_LENGTH)){
		//if(read_pq.empty())
		//prev_mr = merged;
		//else
		  read_pq.push(merged);
	      }
	      else{
		// informative error message!
		if(VERBOSE){
		  cerr << "problem merging read " << read_name << ", splitting read" << endl;
		  cerr << samr.mr << endl;
		  cerr << dangling_mates[read_name].mr << endl;
		  cerr << "To merge, set max segement length (seg_len) higher." << endl;
		}

		// don't throw error for problems merging
		//if(read_pq.empty()){
		// prev_mr = dangling_mates[read_name].mr;
		// read_pq.push(samr.mr);
		//}
		//else{
		  read_pq.push(samr.mr);
		  read_pq.push(dangling_mates[read_name].mr);
		  //}
	      }

	      dangling_mates.erase(read_name);
	    }
	    else{
	      //if(read_pq.empty()){
	      //prev_mr = dangling_mates[read_name].mr;
	      //read_pq.push(samr.mr);
	      //}
	      //else{
		read_pq.push(samr.mr);
		read_pq.push(dangling_mates[read_name].mr);
		//}
	      dangling_mates.erase(read_name);
	    }

	    if(!(read_pq.empty()) &&
	       MappedReadOrderChecker::is_ready(read_pq, samr.mr, MAX_SEGMENT_LENGTH)) {
	      //begin emptying priority queue
	      while(!(read_pq.empty()) &&
		    MappedReadOrderChecker::is_ready(read_pq, samr.mr, MAX_SEGMENT_LENGTH) ){
		empty_pq(prev_mr, read_pq, mapped_reads_file, out);
	      }//end while loop
	    }//end statement for emptying priority queue

	  }
	  else
	    dangling_mates[read_name] = samr;

	}
	else{ 
	    // unmatched, output read
	  if (!samr.is_Trich) revcomp(samr.mr);

	  //if(read_pq.empty())
	  //prev_mr = samr.mr;
	  //else 
	    read_pq.push(samr.mr);

	  if(!(read_pq.empty()) &&
	     MappedReadOrderChecker::is_ready(read_pq, samr.mr, MAX_SEGMENT_LENGTH)) {
	      //begin emptying priority queue
	    while(!(read_pq.empty()) &&
		  MappedReadOrderChecker::is_ready(read_pq, samr.mr, MAX_SEGMENT_LENGTH) ){
	      empty_pq(prev_mr, read_pq, mapped_reads_file, out);
	    }//end while loop
	  }//end statement for emptying priority queue

	}

      // dangling mates is too large, flush dangling_mates of reads
      // on different chroms and too far away 
	if (dangling_mates.size() > MAX_READS_TO_HOLD){
	  
	//   if(VERBOSE)
	//  cerr << "dangling mates too large, emptying" << endl;

	  using std::tr1::unordered_map;
	  unordered_map<string, SAMRecord> tmp;
	  for (unordered_map<string, SAMRecord>::iterator
		 itr = dangling_mates.begin();
	       itr != dangling_mates.end(); ++itr){
	    if (itr->second.mr.r.get_chrom() != samr.mr.r.get_chrom()
		|| (itr->second.mr.r.get_chrom() == samr.mr.r.get_chrom()
		    && itr->second.mr.r.get_end() + MAX_SEGMENT_LENGTH <
		    samr.mr.r.get_start())) {
	      if (!itr->second.is_Trich) revcomp(itr->second.mr);
	      if(itr->second.seg_len >= 0)
		read_pq.push(itr->second.mr);
	    }
	    else
	      tmp[itr->first] = itr->second;
	  }
	  std::swap(tmp, dangling_mates);
	  tmp.clear();

	  if(!(read_pq.empty()) &&
	     MappedReadOrderChecker::is_ready(read_pq, samr.mr, MAX_SEGMENT_LENGTH)) {
	      //begin emptying priority queue
	    while(!(read_pq.empty()) &&
		  MappedReadOrderChecker::is_ready(read_pq, samr.mr, MAX_SEGMENT_LENGTH) ){
	      empty_pq(prev_mr, read_pq, mapped_reads_file, out);
	    }//end while loop
	  }//end statement for emptying priority queue

	}
      
	if (VERBOSE && n_mates % progress_step == 0)
	  cerr << "Processed " << n_mates << " records" << endl;
      }
    }
    // flushing dangling_mates of all remaining ends
    if(!(dangling_mates.empty()) && VERBOSE){
      cerr << "dangling mates not empty" << endl;
      cerr << "dangling_mates.size = " << dangling_mates.size() << endl;
    }
    while (!dangling_mates.empty()){
      if (!dangling_mates.begin()->second.is_Trich)
        revcomp(dangling_mates.begin()->second.mr);
      read_pq.push(dangling_mates.begin()->second.mr);
      dangling_mates.erase(dangling_mates.begin());
    }

    while(!read_pq.empty()){
      empty_pq(prev_mr, read_pq, mapped_reads_file, out);
    }
          
    if (VERBOSE){
      cerr << "Done." << endl;
      cerr << "total mates processed = " << n_mates << endl;
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
