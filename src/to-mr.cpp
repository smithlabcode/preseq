/* to-mr: a program for converting SAM and BAM format to MappedRead
 * format.
 *
 * Copyright (C) 2009-2019 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Meng Zhou, Qiang Song, Andrew Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <cmath>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "htslib_wrapper_deprecated.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::max;
using std::min;
using std::runtime_error;
using std::unordered_map;
using std::swap;


/********Below are functions for merging pair-end reads********/
static void
fill_overlap(const bool pos_str, const MappedRead &mr, const size_t start,
             const size_t end, const size_t offset, string &seq, string &scr) {
  const size_t a = pos_str ? (start - mr.r.get_start()) : (mr.r.get_end() - end);
  const size_t b = pos_str ? (end -  mr.r.get_start()) : (mr.r.get_end() - start);
  copy(mr.seq.begin() + a, mr.seq.begin() + b, seq.begin() + offset);
  copy(mr.scr.begin() + a, mr.scr.begin() + b, scr.begin() + offset);
}

static void
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

  // assert(len > 0);
  // if the above assertion fails, it usually means the mair is
  // discordant (end1 downstream of end2). Also it means the SAM flag
  // of this pair of reads is not properly set. To avoid termination,
  // currently this assertion is ignored but no output will be
  // generated for discordant pairs.

  if (len > 0) {
    string seq(len, 'N');
    string scr(len, 'B');
    if (len <= static_cast<int>(range)) {
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
  }
}

inline static bool
same_read(const size_t suffix_len,
          const MappedRead &a, const MappedRead &b) {
  const string sa(a.r.get_name());
  const string sb(b.r.get_name());
  return std::equal(sa.begin(), sa.end() - suffix_len, sb.begin());
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

static string
get_read_name(const SAMRecord &aln, const size_t suffix_len) {
  return aln.mr.r.get_name().substr(0, aln.mr.r.get_name().size() - suffix_len);
}

int
main(int argc, const char **argv) {

  try {

    string outfile;
    string mapper = "general";
    size_t max_frag_len = 1000;
    size_t max_dangling = 500;
    size_t suffix_len = 1;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Convert the SAM/BAM output from "
                           "to MethPipe mapped read format",
                           "sam/bam_file");
    opt_parse.add_opt("output", 'o', "Name of output file",
                      false, outfile);
    // opt_parse.add_opt("mapper", 'm',
    //                   "Original mapper: bismark, bs_seeker or general",
    //                   true, mapper);
    opt_parse.add_opt("suff", 's', "read name suffix length (default: 1)",
                      false, suffix_len);
    opt_parse.add_opt("max-frag", 'L', "maximum allowed insert size",
                      false, max_frag_len);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc < 3 || opt_parse.help_requested()) {
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
      cerr << "[input file: " << mapped_reads_file << "]" << endl
           << "[output file: "
           << (outfile.empty() ? "stdout" : outfile) << "]" << endl;

    SAMReader_deprecated sam_reader(mapped_reads_file, mapper);
    unordered_map<string, SAMRecord> dangling_mates;

    size_t count = 0;
    SAMRecord aln;

    while (sam_reader >> aln) {
      if (aln.is_mapping_paired) {
        const string read_name(get_read_name(aln, suffix_len));
        if (dangling_mates.find(read_name) != end(dangling_mates)) {
          assert(same_read(suffix_len, aln.mr, dangling_mates[read_name].mr));
          if (aln.is_Trich)
            swap(aln, dangling_mates[read_name]);
          revcomp(aln.mr);

          MappedRead merged;
          int len = 0;
          merge_mates(suffix_len, max_frag_len,
                      dangling_mates[read_name].mr, aln.mr, merged, len);
          if (len <= static_cast<int>(max_frag_len))
            out << merged << endl;
          else if (len > 0)
            out << dangling_mates[read_name].mr << endl << aln.mr << endl;
          dangling_mates.erase(read_name);
        }
        else dangling_mates[read_name] = aln;

        // flush dangling_mates
        if (dangling_mates.size() > max_dangling) {
          unordered_map<string, SAMRecord> tmp;
          for (auto &&mates : dangling_mates)
            if (mates.second.mr.r.get_chrom() < aln.mr.r.get_chrom() ||
                (mates.second.mr.r.get_chrom() == aln.mr.r.get_chrom() &&
                 mates.second.mr.r.get_end() + max_frag_len < aln.mr.r.get_start())) {
              if (!mates.second.is_Trich)
                revcomp(mates.second.mr);
              out << mates.second.mr << endl;
            }
            else tmp[mates.first] = mates.second;
          swap(tmp, dangling_mates);
        }
      }
      else {
        if (!aln.is_Trich)
          revcomp(aln.mr);
        out << aln.mr << endl;
      }
      ++count;
    }
    // flushing dangling_mates
    while (!dangling_mates.empty()) {
      if (!dangling_mates.begin()->second.is_Trich)
        revcomp(dangling_mates.begin()->second.mr);
      out << dangling_mates.begin()->second.mr << endl;
      dangling_mates.erase(dangling_mates.begin());
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
