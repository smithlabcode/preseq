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

#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include <unistd.h>

#include <algorithm>  // std::min
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>  // std::swap
#include <vector>

using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::mt19937;
using std::priority_queue;
using std::runtime_error;
using std::size;
using std::size_t;
using std::string;
using std::uint64_t;
using std::unordered_map;
using std::vector;

//////////////////////////////////////////////////////////////////////
// Data imputation
/////////////////////////////////////////////////////////////////////

static bool
update_pe_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                vector<double> &counts_hist,
                                size_t &current_count) {
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
  else  // next read is same, update current_count
    ++current_count;

  return true;
}

static void
update_se_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                const string &input_file_name,
                                vector<double> &counts_hist,
                                size_t &current_count) {
  // check if reads are sorted
  if (curr_gr.same_chrom(prev_gr) && curr_gr.get_start() < prev_gr.get_start())
    throw runtime_error("locations unsorted in: " + input_file_name);

  if (!curr_gr.same_chrom(prev_gr) ||
      curr_gr.get_start() != prev_gr.get_start()) {
    // next read is new, update counts_hist to include current_count
    // histogram is too small, resize
    if (counts_hist.size() < current_count + 1)
      counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

// comparison function for priority queue

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
  static bool start_check(const GenomicRegion &prev, const GenomicRegion &gr) {
    return (
      chrom_greater(prev, gr) ||
      (prev.same_chrom(gr) && start_greater(prev, gr)) ||
      (prev.same_chrom(gr) && same_start(prev, gr) && end_greater(prev, gr)));
  }
};

typedef priority_queue<GenomicRegion, vector<GenomicRegion>,
                       GenomicRegionOrderChecker>
  ReadPQ;

static bool
is_ready_to_pop(const ReadPQ &pq, const GenomicRegion &gr,
                const size_t max_width) {
  return !pq.top().same_chrom(gr) ||
         pq.top().get_end() + max_width < gr.get_start();
}

static void
empty_pq(GenomicRegion &curr_gr, GenomicRegion &prev_gr, size_t &current_count,
         vector<double> &counts_hist, ReadPQ &read_pq,
         const string &input_file_name) {
  curr_gr = read_pq.top();
  read_pq.pop();

  // update counts hist
  const bool UPDATE_SUCCESS = update_pe_duplicate_counts_hist(
    curr_gr, prev_gr, counts_hist, current_count);
  if (!UPDATE_SUCCESS) {
    std::ostringstream oss;
    oss << "reads unsorted in: " << input_file_name << "\n"
        << "prev = \t" << prev_gr << "\n"
        << "curr = \t" << curr_gr << "\n"
        << "Increase seg_len if in paired end mode";
    throw runtime_error(oss.str());
  }
  prev_gr = curr_gr;
}

/* this code is for BED file input */

size_t
load_counts_BED_se(const string &input_file_name, vector<double> &counts_hist) {
  // resize vals_hist
  counts_hist.clear();
  counts_hist.resize(2, 0.0);

  std::ifstream in(input_file_name);
  if (!in)
    throw runtime_error("problem opening file: " + input_file_name);

  GenomicRegion curr_gr, prev_gr;
  if (!(in >> prev_gr))
    throw runtime_error("problem opening file: " + input_file_name);

  size_t n_reads = 1;
  size_t current_count = 1;
  while (in >> curr_gr) {
    update_se_duplicate_counts_hist(curr_gr, prev_gr, input_file_name,
                                    counts_hist, current_count);
    ++n_reads;
    prev_gr.swap(curr_gr);
  }

  // to account for the last read compared to the one before it.
  if (counts_hist.size() < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return n_reads;
}

size_t
load_counts_BED_pe(const string &input_file_name, vector<double> &counts_hist) {
  // resize vals_hist
  counts_hist.clear();
  counts_hist.resize(2, 0.0);

  std::ifstream in(input_file_name);
  if (!in)
    throw runtime_error("problem opening file: " + input_file_name);

  GenomicRegion curr_gr, prev_gr;
  if (!(in >> prev_gr))
    throw runtime_error("problem opening file: " + input_file_name);

  size_t n_reads = 1;
  size_t current_count = 1;

  // read in file and compare each gr with the one before it
  while (in >> curr_gr) {
    const bool UPDATE_SUCCESS = update_pe_duplicate_counts_hist(
      curr_gr, prev_gr, counts_hist, current_count);
    if (!UPDATE_SUCCESS)
      throw runtime_error("reads unsorted in " + input_file_name);

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
  std::ifstream in(input_file_name);
  if (!in)
    throw runtime_error("problem opening file: " + input_file_name);

  size_t n_counts = 0;
  string buffer;
  while (getline(in, buffer)) {
    if (find(begin(buffer), end(buffer), '\r') != end(buffer))
      throw runtime_error("carriage returns in values file "
                          "(suggests dos or mac formatting)");

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
        n_counts += count;
      }
      else if (val != 0)
        throw runtime_error("problem reading file at line " +
                            toa(n_counts + 1));
    }
    in.peek();
  }
  return n_counts;
}

// returns number of reads from file containing counts histogram
size_t
load_histogram(const string &filename, vector<double> &counts_hist) {
  counts_hist.clear();

  std::ifstream in(filename);
  if (!in)  // if file doesn't open
    throw runtime_error("could not open histogram: " + filename);

  size_t n_reads = 0;
  size_t line_count = 0ul, prev_read_count = 0ul;
  string buffer;
  while (getline(in, buffer)) {
    if (find(begin(buffer), end(buffer), '\r') != end(buffer))
      throw runtime_error("carriage returns in histogram file "
                          "(suggests dos or mac formatting)");

    ++line_count;
    size_t read_count = 0ul;
    double frequency = 0.0;
    std::istringstream is(buffer);
    // error reading input
    if (!(is >> read_count >> frequency))
      throw runtime_error("bad histogram line format:\n" + buffer + "\n" +
                          "(line " + toa(line_count) + ")");

    // histogram is out of order?
    if (read_count < prev_read_count)
      throw runtime_error("bad line order in file " + filename + "\n" +
                          "(line " + toa(line_count) + ")");
    counts_hist.resize(read_count + 1, 0.0);
    counts_hist[read_count] = frequency;
    if (read_count == 0ul) {
      throw runtime_error("counts histograms may not "
                          "include an entry for zero");
    }
    prev_read_count = read_count;
    n_reads += static_cast<size_t>(read_count * frequency);
  }

  return n_reads;
}

/////////////////////////////////////////////////////////
// Loading coverage counts
////////////////////////////////////////////////////////

// probabilistically split genomic regions into mutiple
// genomic regions of width equal to bin_size
static void
SplitGenomicRegion(const GenomicRegion &inputGR, mt19937 &generator,
                   const size_t bin_size, vector<GenomicRegion> &outputGRs) {
  outputGRs.clear();
  GenomicRegion gr(inputGR);

  double frac = static_cast<double>(gr.get_start() % bin_size) / bin_size;
  const size_t width = gr.get_width();

  // ADS: this seems like a bunch of duplicated code just for a single
  // function difference
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  if (dist(generator) > frac) {
    gr.set_start(std::floor(static_cast<double>(gr.get_start()) / bin_size) *
                 bin_size);
    gr.set_end(gr.get_start() + width);
  }
  else {
    gr.set_start(std::ceil(static_cast<double>(gr.get_start()) / bin_size) *
                 bin_size);
    gr.set_end(gr.get_start() + width);
  }

  for (size_t i = 0; i < gr.get_width(); i += bin_size) {
    const size_t curr_start = gr.get_start() + i;
    const size_t curr_end = std::min(gr.get_end(), curr_start + bin_size);
    frac = static_cast<double>(curr_end - curr_start) / bin_size;

    if (dist(generator) <= frac) {
      GenomicRegion binned_gr(gr.get_chrom(), curr_start, curr_start + bin_size,
                              gr.get_name(), gr.get_score(), gr.get_strand());

      outputGRs.push_back(binned_gr);
    }
  }
}

// split a mapped read into multiple genomic regions
// based on the number of bases in each
static void
SplitMappedRead(const MappedRead &inputMR, mt19937 &generator,
                const size_t bin_size, vector<GenomicRegion> &outputGRs) {
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
      const double frac = static_cast<double>(covered_bases) / bin_size;
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      if (dist(generator) <= frac) {
        const size_t curr_start = read_iterator - (read_iterator % bin_size);
        const size_t curr_end = curr_start + bin_size;
        const GenomicRegion binned_gr(
          inputMR.r.get_chrom(), curr_start, curr_end, inputMR.r.get_name(),
          inputMR.r.get_score(), inputMR.r.get_strand());
        outputGRs.push_back(binned_gr);
      }
      total_covered_bases += covered_bases;
      covered_bases = 0;
    }
    seq_iterator++;
    read_iterator++;
  }

  const double frac = static_cast<double>(covered_bases) / bin_size;
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  if (dist(generator) <= frac) {
    const size_t curr_start = read_iterator - (read_iterator % bin_size);
    const size_t curr_end = curr_start + bin_size;
    const GenomicRegion binned_gr(inputMR.r.get_chrom(), curr_start, curr_end,
                                  inputMR.r.get_name(), inputMR.r.get_score(),
                                  inputMR.r.get_strand());
    outputGRs.push_back(binned_gr);
  }
}

size_t
load_coverage_counts_MR(const string &input_file_name, const uint64_t seed,
                        const size_t bin_size, const size_t max_width,
                        vector<double> &coverage_hist) {
  srand(time(0) + getpid());
  // Runif runif(rand());
  std::mt19937 generator(seed);

  std::ifstream in(input_file_name);
  if (!in)
    throw runtime_error("problem opening file: " + input_file_name);

  MappedRead mr;
  if (!(in >> mr))
    throw runtime_error("problem reading from: " + input_file_name);

  // initialize prioirty queue to reorder the split reads
  ReadPQ PQ;

  size_t n_reads = 0;
  size_t n_bins = 0;
  GenomicRegion curr_gr, prev_gr;
  size_t current_count = 1;

  do {
    if (mr.r.get_width() > max_width)
      throw runtime_error("Encountered read of width " + toa(mr.r.get_width()) +
                          "max_width set too small");

    vector<GenomicRegion> splitGRs;
    SplitMappedRead(mr, generator, bin_size, splitGRs);

    n_reads++;
    n_bins += splitGRs.size();

    // add split Genomic Regions to the priority queue
    for (size_t i = 0; i < splitGRs.size(); i++)
      PQ.push(splitGRs[i]);

    // remove Genomic Regions from the priority queue
    if (splitGRs.size() > 0)
      while (!PQ.empty() && is_ready_to_pop(PQ, splitGRs.back(), max_width))
        empty_pq(curr_gr, prev_gr, current_count, coverage_hist, PQ,
                 input_file_name);
  } while (in >> mr);

  // done adding reads, now spit the rest out
  while (!PQ.empty())
    empty_pq(curr_gr, prev_gr, current_count, coverage_hist, PQ,
             input_file_name);

  return n_reads;
}

size_t
load_coverage_counts_GR(const string &input_file_name, const uint64_t seed,
                        const size_t bin_size, const size_t max_width,
                        vector<double> &coverage_hist) {
  srand(time(0) + getpid());
  std::mt19937 generator(seed);

  std::ifstream in(input_file_name);
  if (!in)
    throw runtime_error("problem opening file: " + input_file_name);

  GenomicRegion inputGR;
  if (!(in >> inputGR))
    throw runtime_error("problem reading from: " + input_file_name);

  // initialize prioirty queue to reorder the split reads
  ReadPQ PQ;

  // prev and current Genomic Regions to compare
  GenomicRegion curr_gr, prev_gr;
  size_t n_reads = 0;
  size_t current_count = 1;

  do {
    vector<GenomicRegion> splitGRs;
    SplitGenomicRegion(inputGR, generator, bin_size, splitGRs);

    // add split Genomic Regions to the priority queue
    for (size_t i = 0; i < splitGRs.size(); i++)
      PQ.push(splitGRs[i]);

    if (splitGRs.size() > 0) {
      // remove Genomic Regions from the priority queue
      while (!PQ.empty() && is_ready_to_pop(PQ, splitGRs.back(), max_width))
        empty_pq(curr_gr, prev_gr, current_count, coverage_hist, PQ,
                 input_file_name);
    }
    n_reads++;
  } while (in >> inputGR);

  // done adding reads, now spit the rest out
  while (!PQ.empty())
    empty_pq(curr_gr, prev_gr, current_count, coverage_hist, PQ,
             input_file_name);

  return n_reads;
}

#ifdef HAVE_HTSLIB
/* Deal with SAM/BAM format only if we have htslib
 */

#include "bam_record_utils.hpp"  // from dnmtools

#include <bamxx.hpp>  // from bamxx
#include <htslib_wrapper.hpp>

static inline void
swap(bamxx::bam_rec &a, bamxx::bam_rec &b) {
  std::swap(a.b, b.b);
}

static inline void
update_se_duplicate_counts_hist_BAM(
  const int32_t curr_tid, const hts_pos_t curr_pos, const int32_t prev_tid,
  const hts_pos_t prev_pos, const string &inputfile,
  vector<double> &counts_hist, size_t &current_count) {
  // check if reads are sorted; situation of unordered chroms is taken
  // care of outside this function
  if (curr_tid == prev_tid && curr_pos < prev_pos)
    throw runtime_error("locations unsorted in: " + inputfile);

  if (curr_tid != prev_tid || curr_pos != prev_pos) {
    // next read is new, update counts_hist to include current_count
    if (counts_hist.size() < current_count + 1) {
      // histogram is too small, resize
      counts_hist.resize(current_count + 1, 0.0);
    }
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

size_t
load_counts_BAM_se(/*const size_t n_threads, */
                   const string &inputfile, vector<double> &counts_hist) {
  // bamxx::bam_tpool tp(n_threads);

  bamxx::bam_in hts(inputfile);  // assume already checked
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw runtime_error("failed to read header");

  // if (n_threads > 1) tp.set_io(hts);

  // find first mapped read to start
  bamxx::bam_rec aln;
  while (hts.read(hdr, aln) && get_tid(aln) == -1)
    ;

  size_t n_reads{};
  // if all reads unmapped, must return
  if (get_tid(aln) == -1)
    return n_reads;

  // to check that reads are sorted properly
  vector<bool> chroms_seen(get_n_targets(hdr), false);

  // start with prev_aln being first read
  bamxx::bam_rec prev_aln;
  swap(aln, prev_aln);
  int32_t prev_tid = get_tid(prev_aln);
  int32_t prev_pos = get_pos(prev_aln);

  // start with count of 1 for first read seen
  size_t current_count = 1;

  while (hts.read(hdr, aln)) {
    const int32_t curr_tid = get_tid(aln);
    const int32_t curr_pos = get_pos(aln);
    if (curr_tid != prev_tid) {
      if (chroms_seen[curr_tid])
        throw runtime_error("input not sorted");
      chroms_seen[curr_tid] = true;
    }

    // check that mapped read is not secondary
    if (curr_tid != -1) {  // check that read is mapped
      update_se_duplicate_counts_hist_BAM(curr_tid, curr_pos, prev_tid,
                                          prev_pos, inputfile, counts_hist,
                                          current_count);
      ++n_reads;
    }
    prev_pos = curr_pos;
    prev_tid = curr_tid;
  }

  // to account for the last read compared to the one before it.
  if (size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return n_reads;
}

static inline void
update_pe_duplicate_counts_hist_BAM(
  // clang-format off
  const int32_t curr_tid, const hts_pos_t curr_pos,
  const int32_t curr_mtid, const hts_pos_t curr_mpos,  // vals for curr
  const int32_t prev_tid, const hts_pos_t prev_pos,
  const int32_t prev_mtid, const hts_pos_t prev_mpos,  // vals for prev
  // clang-format on
  const string &inputfile, vector<double> &counts_hist, size_t &current_count) {
  // check if reads are sorted; situation of unordered chroms is taken
  // care of outside this function
  if (curr_tid == prev_tid && curr_pos < prev_pos)
    throw runtime_error("locations unsorted in: " + inputfile);

  if (curr_tid != prev_tid || curr_pos != prev_pos || curr_mtid != prev_mtid ||
      curr_mpos != prev_mpos) {
    // next read is new, update counts_hist to include current_count
    if (counts_hist.size() < current_count + 1) {
      // histogram is too small, resize
      counts_hist.resize(current_count + 1, 0.0);
    }
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

size_t
load_counts_BAM_pe(const string &inputfile, vector<double> &counts_hist) {
  bamxx::bam_in hts(inputfile);  // assume already checked
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw runtime_error("failed to read header");

  // if (n_threads > 1) tp.set_io(hts);

  // find first mapped read to start
  bamxx::bam_rec aln;
  while (hts.read(hdr, aln) && get_tid(aln) == -1)
    ;

  size_t n_reads{};
  // if all reads unmapped, must return
  if (get_tid(aln) == -1)
    return n_reads;

  // to check that reads are sorted properly
  vector<bool> chroms_seen(get_n_targets(hdr), false);

  // start with prev_aln being first read
  bamxx::bam_rec prev_aln;
  swap(aln, prev_aln);

  int32_t prev_tid = get_tid(prev_aln);
  int32_t prev_pos = get_pos(prev_aln);
  int32_t prev_mtid = get_mtid(prev_aln);
  int32_t prev_mpos = get_mpos(prev_aln);

  // start with count of 1 for first read seen
  size_t current_count = 1;

  while (hts.read(hdr, aln)) {
    const int32_t curr_tid = get_tid(aln);
    const int32_t curr_pos = get_pos(aln);
    const int32_t curr_mtid = get_mtid(aln);
    const int32_t curr_mpos = get_mpos(aln);

    if (curr_mtid < curr_tid ||
        (curr_mtid == curr_tid && curr_mpos == curr_mtid)) {
      // ADS: this means that we should already have processed this
      // same read pair earlier
      continue;
    }

    if (curr_tid != prev_tid) {
      if (chroms_seen[curr_tid])
        throw runtime_error("input not sorted");
      chroms_seen[curr_tid] = true;
    }

    // check that mapped read is not secondary
    if (curr_tid != -1) {  // check that read is mapped
      update_pe_duplicate_counts_hist_BAM(
        curr_tid, curr_pos, curr_mtid, curr_mpos,  // vals for curr
        prev_tid, prev_pos, prev_mtid, prev_mpos,  // vals for prev
        inputfile, counts_hist, current_count);

      ++n_reads;
    }
    prev_tid = curr_tid;
    prev_pos = curr_pos;
    prev_mtid = curr_mtid;
    prev_mpos = curr_mpos;
  }

  // to account for the last read compared to the one before it.
  if (size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return n_reads;
}

#endif
