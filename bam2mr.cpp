/*    bam2mr: program to BAM format to smithlab's MappedRead format
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith and
 *                       Timothy Daley
 *
 *    Authors: Andrew D. Smith & Timothy Daley 
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <sys/types.h>
#include <unistd.h>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MappedRead.hpp"

#include "api/BamReader.h"
#include "api/BamAlignment.h"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

using std::tr1::unordered_map;

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;

static void
BamAlignmentToMappedRead(const unordered_map<size_t, string> &chrom_lookup,
			 const BamAlignment &ba,
			 MappedRead &mr) {

  const unordered_map<size_t, string>::const_iterator 
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));
  
  const string chrom = the_chrom->second;
  size_t start, end;
  if(ba.Position + ba.InsertSize < 0){
    start = ba.Position;
    end = ba.Position - ba.InsertSize;
  }
  else{
    start = std::min(ba.Position, ba.Position + ba.InsertSize);
    end = std::max(ba.Position, ba.Position + ba.InsertSize);
  }
  const string name(ba.Name);
  const float score = ba.MapQuality;
  const char strand = (ba.IsReverseStrand() ? '-' : '+');
  const string seq = ba.AlignedBases;
  const string scr = ba.Qualities;
  mr.r = GenomicRegion(chrom, start, end, name, score, strand);
  mr.seq = seq;
  mr.scr = scr;
}



static void
fill_overlap(const bool pos_str, const MappedRead &mr, const size_t start, 
	     const size_t end, const size_t offset, string &seq, string &scr) {
  const size_t a = pos_str ? (start - mr.r.get_start()) : (mr.r.get_end() - end);
  const size_t b = pos_str ? (end -  mr.r.get_start()) : (mr.r.get_end() - start);
  copy(mr.seq.begin() + a, mr.seq.begin() + b, seq.begin() + offset);
  copy(mr.scr.begin() + a, mr.scr.begin() + b, scr.begin() + offset);
}



static bool
merge_mates(const bool VERBOSE, const size_t range, const MappedRead &one, 
	    const MappedRead &two, MappedRead &merged) {
  
  const bool pos_str = one.r.pos_strand();
  const size_t overlap_start = std::max(one.r.get_start(), two.r.get_start());
  const size_t overlap_end = std::min(one.r.get_end(), two.r.get_end());

  const size_t one_left = pos_str ? 
    one.r.get_start() : std::max(overlap_end, one.r.get_start());
  const size_t one_right = 
    pos_str ? std::min(overlap_start, one.r.get_end()) : one.r.get_end();
  
  const size_t two_left = pos_str ? 
    std::max(overlap_end, two.r.get_start()) : two.r.get_start();
  const size_t two_right = pos_str ? 
    two.r.get_end() : std::min(overlap_start, two.r.get_end());
  
  const int len = pos_str ? (two_right - one_left) : (one_right - two_left);
 
   
  if(len < 0){
    if(VERBOSE){
      cerr << one << endl;
      cerr << two << endl;
      cerr << "len = " << len << endl;
    } 
    return false;
  }
  
  if(!(one_left <= one_right && two_left <= two_right)){
    cerr << one << endl;
    cerr << two << endl;
  }

  assert(one_left <= one_right && two_left <= two_right);
  assert(overlap_start >= overlap_end || static_cast<size_t>(len) == 
	 ((one_right - one_left) + (two_right - two_left) + (overlap_end - overlap_start)));
  
  string seq(len, 'N');
  string scr(len, 'B');
  if (len > 0 && len <= static_cast<int>(range)) {
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
  merged.r.set_name("FRAG:" + name.substr(0, name.size()));

  return true;
}

static bool
BamAlignmentPeToMappedRead(const bool VERBOSE,
			   const unordered_map<size_t, string> &chrom_lookup,
			   const BamAlignment &ba_1,
			   const BamAlignment &ba_2, 
			   const size_t max_frag_len,
			   MappedRead &mr) {

  MappedRead mr_1, mr_2, merged_mr;
  BamAlignmentToMappedRead(chrom_lookup, ba_1, mr_1);
  BamAlignmentToMappedRead(chrom_lookup, ba_2, mr_2);

  return merge_mates(VERBOSE, max_frag_len, mr_1, mr_2, mr);

}

inline static bool
same_read(const BamAlignment &a, const BamAlignment &b) {
  const string sa(a.Name);
  const string sb(b.Name);
  return std::equal(sa.begin(), sa.end(), sb.begin());
}



int 
main(int argc, const char **argv)  {
  try {

    bool VERBOSE = false;
    string outfile;
    size_t max_frag_len = 10000;
      
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "convert paired end BAM format reads to mapped reads format");
    opt_parse.add_opt("fraglen", 'L', "max fragment length", false, max_frag_len);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("outfile", 'o', "output file name", false, outfile);
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
    if (leftover_args.size() != 1) {
      cerr << "input file missing" << endl;
      return EXIT_SUCCESS;
    }
    const string infile(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    BamReader reader;
    reader.Open(infile);
    
    // Get header and reference
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    unordered_map<size_t, string> chrom_lookup;
    for (size_t i = 0; i < refs.size(); ++i)
      chrom_lookup[i] = refs[i].RefName;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    BamAlignment bam_1, bam_2;
    while (reader.GetNextAlignment(bam_1)) {
      MappedRead mr;

      // if the reads are a proper pair, get next alignment (should be the mate)
      // and merge alignment
      if(bam_1.IsPaired() && bam_1.IsMapped() && 
	 bam_1.IsProperPair() && bam_1.IsPrimaryAlignment()){
	do{
	  reader.GetNextAlignment(bam_2);
	} while (!(bam_2.IsPaired() && bam_2.IsMapped() && 
		   bam_2.IsProperPair() && bam_2.IsPrimaryAlignment()));
	// if the next read is not the mate, they will not have the same name
	if(!same_read(bam_1, bam_2)){
	  MappedRead mr_1, mr_2;
	  BamAlignmentToMappedRead(chrom_lookup, bam_1, mr_1);
	  BamAlignmentToMappedRead(chrom_lookup, bam_2, mr_2);
	  const string sa(bam_1.Name);
	  const string sb(bam_2.Name);
	  cerr << sa << '\t' << mr_1 << endl;
	  cerr << sb << '\t' << mr_2 << endl;
	  throw SMITHLABException("Reads not sorted by name");
	}
	bool merge_success = BamAlignmentPeToMappedRead(VERBOSE, chrom_lookup, bam_1, bam_2, max_frag_len, mr);
	if(merge_success)
	  out << mr << endl;
	/*
	else{
	    MappedRead mr_1, mr_2;
	    BamAlignmentToMappedRead(chrom_lookup, bam_1, mr_1);
	    BamAlignmentToMappedRead(chrom_lookup, bam_1, mr_2);
	    cerr << mr_1 << endl;
	    cerr << mr_2 << endl;
	    throw SMITHLABException("Problem merging mates");
	  }
	*/
      }
      // if the read is not a proper pair, convert it mapped reads
      else if(bam_1.IsMapped() && bam_1.IsPrimaryAlignment()){
	BamAlignmentToMappedRead(chrom_lookup, bam_1, mr);
	out << mr << endl;
      }

      // if the read is not mapped or is not primary alignment, do nothing

    }
    reader.Close();
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
