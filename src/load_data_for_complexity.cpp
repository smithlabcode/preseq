/* Copyright (C) 2013-2026 Andrew D. Smith and Timothy Daley
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "load_data_for_complexity.hpp"

#include "bamxx/bamxx.hpp"

#include <htslib/sam.h>

#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include <unistd.h>

#include <algorithm>  // std::min
#include <cassert>
#include <iostream>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_map>
#include <utility>  // std::swap

[[nodiscard]] static auto
update_pe_duplicate_counts_hist(const auto &curr, const auto &prev,
                                std::vector<double> &counts_hist,
                                std::size_t &current_count) -> bool {
  // check if reads are sorted
  if (curr.same_chrom(prev) && curr.get_start() < prev.get_start() &&
      curr.get_end() < prev.get_end()) {
    return false;
  }

  // check if next read is new, and if so update counts_hist to
  // include current_count
  if (!curr.same_chrom(prev) || curr.get_start() != prev.get_start() ||
      curr.get_end() != prev.get_end()) {
    // histogram is too small, resize
    if (std::size(counts_hist) < current_count + 1)
      counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;

  return true;
}

static void
update_se_duplicate_counts_hist(const auto &curr, const auto &prev,
                                const std::string input_file_name,
                                std::vector<double> &counts_hist,
                                std::size_t &current_count) {
  // check if reads are sorted
  if (curr.same_chrom(prev) && curr.get_start() < prev.get_start())
    throw std::runtime_error("locations unsorted in: " + input_file_name);

  if (!curr.same_chrom(prev) || curr.get_start() != prev.get_start())
  // next read is new, update counts_hist to include current_count
  {
    // histogram is too small, resize
    if (std::size(counts_hist) < current_count + 1)
      counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

// comparison function for priority queue

/**************** FOR CLARITY BELOW WHEN COMPARING READS *************/
[[nodiscard]] static inline auto
chrom_greater(const auto &a, const auto &b) -> bool {
  return a.get_chrom() > b.get_chrom();
}
[[nodiscard]] static inline auto
same_start(const auto &a, const auto &b) -> bool {
  return a.get_start() == b.get_start();
}
[[nodiscard]] static inline auto
start_greater(const auto &a, const auto &b) -> bool {
  return a.get_start() > b.get_start();
}
[[nodiscard]] static inline auto
end_greater(const auto &a, const auto &b) -> bool {
  return a.get_end() > b.get_end();
}
/******************************************************************************/

struct GenomicRegionOrderChecker {
  [[nodiscard]] auto operator()(const GenomicRegion &prev,
                                const GenomicRegion &gr) const -> bool {
    return start_check(prev, gr);
  }
  [[nodiscard]] static auto start_check(const GenomicRegion &prev,
                                        const GenomicRegion &gr) -> bool {
    return chrom_greater(prev, gr) ||                           //
           (prev.same_chrom(gr) && start_greater(prev, gr)) ||  //
           (prev.same_chrom(gr) && same_start(prev, gr) &&
            end_greater(prev, gr));
  }
};

using ReadPQ = std::priority_queue<GenomicRegion, std::vector<GenomicRegion>,
                                   GenomicRegionOrderChecker>;

[[nodiscard]] static auto
is_ready_to_pop(const ReadPQ &pq, const auto &gr,
                const std::size_t max_width) -> bool {
  return !pq.top().same_chrom(gr) ||
         pq.top().get_end() + max_width < gr.get_start();
}

static void
empty_pq(auto &curr, auto &prev, std::size_t &current_count,
         std::vector<double> &counts_hist, ReadPQ &read_pq,
         const std::string &input_file_name) {
  curr = read_pq.top();
  read_pq.pop();

  // update counts hist
  const bool update_success =
    update_pe_duplicate_counts_hist(curr, prev, counts_hist, current_count);
  if (!update_success) {
    std::ostringstream oss;
    oss << "reads unsorted in: " << input_file_name << "\n"
        << "prev = \t" << prev << "\n"
        << "curr = \t" << curr << "\n"
        << "Increase seg_len if in paired end mode";
    throw std::runtime_error(oss.str());
  }
  prev = curr;
}

/*
 * This code is used to deal with read data in BAM format.
 */
#ifdef HAVE_HTSLIB
// switching dependency on bamtools to samtools
// #include "htslib_wrapper_deprecated.hpp"

[[nodiscard]] auto
load_counts_BAM_se(const std::string &input_file_name,
                   std::vector<double> &counts_hist) -> std::size_t {
  // const std::string mapper = "general";
  // SAMReader_deprecated sam_reader(input_file_name, mapper);
  // if (!sam_reader)
  //   throw std::runtime_error("problem opening input file " +
  //   input_file_name);

  // open the hts SAM/BAM input file and get the header
  bamxx::bam_in hts(input_file_name);
  if (!hts)
    throw std::runtime_error("failed to open input file");
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  bamxx::bam_rec aln;

  // SAMRecord samr;
  // sam_reader >> samr;
  std::size_t n_reads = 1;
  // resize vals_hist, make sure it starts out empty
  counts_hist.clear();
  counts_hist.resize(2, 0.0);
  std::size_t current_count = 1;

  const auto get_seq = [](const auto &aln) {
    auto qlen = aln.b->core.l_qseq;
    auto seq = bam_get_seq(aln.b);
    std::string seq_str;
    for (auto i = 0; i < qlen; ++i)
      seq_str += seq_nt16_str[bam_seqi(seq, i)];
    return seq_str;
  };

  const auto get_mr = [&](const auto &aln) {
    MappedRead mr;
    // clang-format off
    mr.r = GenomicRegion(
      sam_hdr_tid2name(hdr.h, aln.b->core.tid),
      aln.b->core.pos,
      bam_endpos(aln.b),
      bam_get_qname(aln.b),
      0.0,
      bam_is_rev(aln.b) ? '-' : '+'
    );
    // clang-format on
    mr.seq = get_seq(aln);
    mr.scr.resize(aln.b->core.l_qseq, 'B');
    return mr;
  };
  //   return MappedRead{
  //     GenomicRegion(sam_hdr_tid2name(hdr.h, aln.b->core.tid),
  //                     aln.b->core.pos, bam_endpos(aln.b),
  //                     bam_get_qname(aln.b), 0.0,
  //                     bam_is_rev(aln.b) ? '-' : '+'),
  //       get_seq(aln),
  //       std::string(aln.b->core.l_qseq, 'B')
  //     };
  //   // clang-format on
  // };

  //   /*! @abstract the read is paired in sequencing, no matter whether it is
  //   mapped in a pair */
  // #define BAM_FPAIRED        1
  // /*! @abstract the read is mapped in a proper pair */
  // #define BAM_FPROPER_PAIR   2
  // /*! @abstract the read itself is unmapped; conflictive with
  // BAM_FPROPER_PAIR */ #define BAM_FUNMAP         4
  // /*! @abstract the mate is unmapped */
  // #define BAM_FMUNMAP        8
  // /*! @abstract the read is mapped to the reverse strand */
  // #define BAM_FREVERSE      16
  // /*! @abstract the mate is mapped to the reverse strand */
  // #define BAM_FMREVERSE     32
  // /*! @abstract this is read1 */
  // #define BAM_FREAD1        64
  // /*! @abstract this is read2 */
  // #define BAM_FREAD2       128
  // /*! @abstract not primary alignment */
  // #define BAM_FSECONDARY   256
  // /*! @abstract QC failure */
  // #define BAM_FQCFAIL      512
  // /*! @abstract optical or PCR duplicate */
  // #define BAM_FDUP        1024
  // /*! @abstract supplementary alignment */
  // #define BAM_FSUPPLEMENTARY 2048

  constexpr auto check_flag = [&](const auto &aln, const auto flag) {
    return (aln.b->core.flag & flag) != 0;
  };
  constexpr auto get_tid = [](const auto &aln) { return aln.b->core.tid; };

  MappedRead prev_mr, curr_mr;
  prev_mr = get_mr(aln);  // samr.mr;

  while (hts.read(hdr, aln)) {
    const std::int32_t tid = get_tid(aln);
    if (tid == -1)  // ADS: skip unmapped lines reads
      continue;
    // while (sam_reader >> samr) {
    // only convert mapped and primary reads
    if (!check_flag(aln, BAM_FUNMAP) && !check_flag(aln, BAM_FSECONDARY) &&
        (!check_flag(aln, BAM_FPROPER_PAIR) || check_flag(aln, BAM_FREAD1))) {
      // if (samr.is_primary && samr.is_mapped) {
      //   // ignore unmapped reads & secondary alignments
      //   if (!samr.is_mapping_paired || samr.is_mapping_paired &&
      //   samr.is_Trich) {
      // only count unpaired reads or the left mate of paired reads
      curr_mr = get_mr(aln);  // samr.mr;
      update_se_duplicate_counts_hist(curr_mr.r, prev_mr.r, input_file_name,
                                      counts_hist, current_count);
      // update number of reads and prev read
      ++n_reads;
      prev_mr = curr_mr;  // samr.mr;
      // }
    }
  }

  // to account for the last read compared to the one before it.
  if (std::size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return n_reads;
}

/********Below are functions for merging pair-end reads********/

static bool
merge_mates(const std::size_t suffix_len, const auto &one, const auto &two,
            auto &merged, int &len) {
  assert(one.same_chrom(two));
  const std::size_t read_start = std::min(one.get_start(), two.get_start());
  const std::size_t read_end = std::max(one.get_end(), two.get_end());

  len = read_end - read_start;

  if (len < 0)
    return false;

  merged = one;
  merged.set_start(read_start);
  merged.set_end(read_end);
  merged.set_score(one.get_score() + two.get_score());

  const std::string name(one.get_name());
  merged.set_name("FRAG:" + name.substr(0, std::size(name) - suffix_len));

  return true;
}

[[nodiscard]] static inline auto
same_read(const std::size_t suffix_len, const MappedRead &a,
          const MappedRead &b) -> bool {
  const auto sa = a.r.get_name();
  const auto sb = b.r.get_name();
  bool SAME_NAME = false;
  if (sa == sb)
    SAME_NAME = true;
  return (SAME_NAME && a.r.same_chrom(b.r));
}

// return true if the genomic region is null
[[nodiscard]] static inline auto
GenomicRegionIsNull(const GenomicRegion &gr) -> bool {
  GenomicRegion null_gr;
  return gr == null_gr
}

static void
empty_pq(GenomicRegion &prev,
         std::priority_queue<GenomicRegion, std::vector<GenomicRegion>,
                             GenomicRegionOrderChecker> &read_pq,
         const std::string &input_file_name, std::vector<double> &counts_hist,
         std::size_t &current_count) {
  GenomicRegion curr = read_pq.top();
  read_pq.pop();

  // check if reads are sorted
  if (curr.same_chrom(prev) && curr.get_start() < prev.get_start() &&
      curr.get_end() < prev.get_end()) {
    std::ostringstream oss;
    oss << "reads unsorted in: " << input_file_name << "\n"
        << "prev = \t" << prev << "\n"
        << "curr = \t" << curr << "\n"
        << "Increase seg_len if in paired end mode";
    throw std::runtime_error(oss.str());
  }

  if (GenomicRegionIsNull(prev))
    current_count = 1;
  else {
    std::ostringstream oss;
    bool UPDATE_HIST =
      update_pe_duplicate_counts_hist(curr, prev, counts_hist, current_count);
    if (!UPDATE_HIST) {
      oss << "locations unsorted in: " << input_file_name << "\n"
          << "prev = \t" << prev << "\n"
          << "curr = \t" << curr << "\n";
      throw std::runtime_error(oss.str());
    }
  }
  prev = curr;
}

[[nodiscard]] auto
load_counts_BAM_pe(const std::string &input_file_name,
                   const std::size_t MAX_SEGMENT_LENGTH,
                   const std::size_t MAX_READS_TO_HOLD, std::size_t &n_paired,
                   std::size_t &n_mates,
                   std::vector<double> &counts_hist) -> std::size_t {
  const std::string mapper = "general";
  SAMReader_deprecated sam_reader(input_file_name, mapper);

  // check sam_reader
  if (!sam_reader)
    throw std::runtime_error("problem opening input file " + input_file_name);

  SAMRecord samr;
  // resize vals_hist, make sure it starts out empty
  counts_hist.clear();
  counts_hist.resize(2, 0.0);
  std::size_t current_count = 0;
  std::size_t suffix_len = 0;
  n_paired = 0;
  n_mates = 0;
  std::size_t n_unpaired = 0;
  std::size_t progress_step = 1000000;

  GenomicRegion prev;

  std::priority_queue<GenomicRegion, std::vector<GenomicRegion>,
                      GenomicRegionOrderChecker>
    read_pq;

  std::unordered_map<std::string, SAMRecord> dangling_mates;

  while (sam_reader >> samr) {
    // only convert mapped and primary reads
    if (samr.is_primary && samr.is_mapped) {
      ++n_mates;

      // deal with paired-end stuff
      if (samr.is_mapping_paired) {
        const std::size_t name_len = samr.mr.r.get_name().size() - suffix_len;
        const std::string read_name(samr.mr.r.get_name().substr(0, name_len));

        if (dangling_mates.find(read_name) != std::cend(dangling_mates)) {
          // other end is in dangling mates, merge the two mates
          if (same_read(suffix_len, samr.mr, dangling_mates[read_name].mr)) {
            if (samr.is_Trich)
              std::swap(samr, dangling_mates[read_name]);
            GenomicRegion merged;
            int len = 0;
            const bool MERGE_SUCCESS = merge_mates(
              suffix_len, MAX_SEGMENT_LENGTH, dangling_mates[read_name].mr.r,
              samr.mr.r, merged, len);
            // merge success!
            if (MERGE_SUCCESS && len >= 0 &&
                len <= static_cast<int>(MAX_SEGMENT_LENGTH)) {
              read_pq.push(merged);
              ++n_paired;
            }
            else {
              // informative error message!
              // if (VERBOSE) {
              std::cerr << "problem merging read " << read_name
                        << ", splitting read" << std::endl
                        << samr.mr << "\t" << samr.is_mapping_paired
                        << std::endl
                        << dangling_mates[read_name].mr << "\t"
                        << dangling_mates[read_name].is_mapping_paired
                        << std::endl
                        << "To merge, set max segement "
                        << "length (seg_len) higher." << std::endl;
              // }
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
        else  // didn't find read in dangling_mates, store for later
          dangling_mates[read_name] = samr;
      }
      else {
        read_pq.push(samr.mr.r);
        ++n_unpaired;
      }

      // dangling mates is too large, flush dangling_mates of reads
      // on different chroms and too far away
      if (std::size(dangling_mates) > MAX_READS_TO_HOLD) {
        std::unordered_map<std::string, SAMRecord> tmp;
        for (auto itr = std::begin(dangling_mates);
             itr != std::end(dangling_mates); ++itr) {
          if (itr->second.mr.r.get_chrom() != samr.mr.r.get_chrom() ||
              (itr->second.mr.r.get_chrom() == samr.mr.r.get_chrom() &&
               itr->second.mr.r.get_end() + MAX_SEGMENT_LENGTH <
                 samr.mr.r.get_start())) {
            if (itr->second.seg_len >= 0) {
              read_pq.push(itr->second.mr.r);
              ++n_unpaired;
            }
          }
          else
            tmp[itr->first] = itr->second;
        }
        std::swap(tmp, dangling_mates);
        tmp.clear();
      }

      // now empty the priority queue
      if (!(read_pq.empty()) &&
          is_ready_to_pop(read_pq, samr.mr.r, MAX_SEGMENT_LENGTH)) {
        // begin emptying priority queue
        while (!(read_pq.empty()) &&
               is_ready_to_pop(read_pq, samr.mr.r, MAX_SEGMENT_LENGTH)) {
          empty_pq(prev, read_pq, input_file_name, counts_hist, current_count);
        }
      }

      if (/*VERBOSE && */ n_mates % progress_step == 0)
        std::cerr << "Processed " << n_mates << " records" << std::endl;
    }
  }

  // empty dangling mates of any excess reads
  while (!dangling_mates.empty()) {
    read_pq.push(std::begin(dangling_mates)->second.mr.r);
    dangling_mates.erase(std::begin(dangling_mates));
    ++n_unpaired;
  }

  // final iteration
  while (!read_pq.empty())
    empty_pq(prev, read_pq, input_file_name, counts_hist, current_count);

  if (std::size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);

  ++counts_hist[current_count];

  assert(read_pq.empty());

  const std::size_t n_reads = n_unpaired + n_paired;

  // if (VERBOSE)
  std::cerr << "paired = " << n_paired << std::endl
            << "unpaired = " << n_unpaired << std::endl;

  return n_reads;
}

#endif

/* this code is for BED file input */

auto
load_counts_BED_se(const std::string input_file_name,
                   std::vector<double> &counts_hist) -> std::size_t {
  // resize vals_hist
  counts_hist.clear();
  counts_hist.resize(2, 0.0);

  std::ifstream in(input_file_name);
  if (!in)
    throw std::runtime_error("problem opening file: " + input_file_name);

  GenomicRegion curr, prev;
  if (!(in >> prev))
    throw std::runtime_error("problem opening file: " + input_file_name);

  std::size_t n_reads = 1;
  std::size_t current_count = 1;
  while (in >> curr) {
    update_se_duplicate_counts_hist(curr, prev, input_file_name, counts_hist,
                                    current_count);
    ++n_reads;
    prev.swap(curr);
  }

  // to account for the last read compared to the one before it.
  if (std::size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return n_reads;
}

auto
load_counts_BED_pe(const std::string input_file_name,
                   std::vector<double> &counts_hist) -> std::size_t {
  // resize vals_hist
  counts_hist.clear();
  counts_hist.resize(2, 0.0);

  std::ifstream in(input_file_name);
  if (!in)
    throw std::runtime_error("problem opening file: " + input_file_name);

  GenomicRegion curr, prev;
  if (!(in >> prev))
    throw std::runtime_error("problem opening file: " + input_file_name);

  std::size_t n_reads = 1;
  std::size_t current_count = 1;

  // read in file and compare each gr with the one before it
  while (in >> curr) {
    const bool UPDATE_SUCCESS =
      update_pe_duplicate_counts_hist(curr, prev, counts_hist, current_count);
    if (!UPDATE_SUCCESS)
      throw std::runtime_error("reads unsorted in " + input_file_name);

    ++n_reads;
    prev.swap(curr);
  }

  if (std::size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);

  // to account for the last read compared to the one before it.
  ++counts_hist[current_count];

  return n_reads;
}

/* text file input */
auto
load_counts(const std::string &input_file_name,
            std::vector<double> &counts_hist) -> std::size_t {
  std::ifstream in(input_file_name);
  if (!in)
    throw std::runtime_error("problem opening file: " + input_file_name);

  std::size_t n_counts = 0;
  std::string buffer;
  while (getline(in, buffer)) {
    if (find(begin(buffer), end(buffer), '\r') != end(buffer))
      throw std::runtime_error("carriage returns in values file "
                               "(suggests dos or mac formatting)");

    std::istringstream iss(buffer);
    if (iss.good()) {
      double val;
      iss >> val;
      if (val > 0) {
        const auto count = static_cast<std::size_t>(val);
        // histogram is too small, resize
        if (std::size(counts_hist) < count + 1)
          counts_hist.resize(count + 1, 0.0);
        ++counts_hist[count];
        n_counts += count;
      }
      else if (val != 0)
        throw std::runtime_error("problem reading file at line " +
                                 toa(n_counts + 1));
    }
    in.peek();
  }
  return n_counts;
}

// returns number of reads from file containing counts histogram
auto
load_histogram(const std::string &filename,
               std::vector<double> &counts_hist) -> std::size_t {
  counts_hist.clear();

  std::ifstream in(filename);
  if (!in)  // if file doesn't open
    throw std::runtime_error("could not open histogram: " + filename);

  std::size_t n_reads = 0;
  std::size_t line_count = 0ul, prev_read_count = 0ul;
  std::string buffer;
  while (getline(in, buffer)) {
    if (find(begin(buffer), end(buffer), '\r') != end(buffer))
      throw std::runtime_error("carriage returns in histogram file "
                               "(suggests dos or mac formatting)");

    ++line_count;
    std::size_t read_count = 0ul;
    double frequency = 0.0;
    std::istringstream is(buffer);
    // error reading input
    if (!(is >> read_count >> frequency))
      throw std::runtime_error("bad histogram line format:\n" + buffer + "\n" +
                               "(line " + toa(line_count) + ")");

    // histogram is out of order?
    if (read_count < prev_read_count)
      throw std::runtime_error("bad line order in file " + filename + "\n" +
                               "(line " + toa(line_count) + ")");
    counts_hist.resize(read_count + 1, 0.0);
    counts_hist[read_count] = frequency;
    if (read_count == 0ul) {
      throw std::runtime_error("counts histograms may not "
                               "include an entry for zero");
    }
    prev_read_count = read_count;
    n_reads += static_cast<std::size_t>(read_count * frequency);
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
  const std::size_t width = gr.get_width();

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

  for (std::size_t i = 0; i < gr.get_width(); i += bin_size) {
    const std::size_t curr_start = gr.get_start() + i;
    const std::size_t curr_end = std::min(gr.get_end(), curr_start + bin_size);
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
SplitMappedRead(const MappedRead &inputMR, std::mt19937 &generator,
                const std::size_t bin_size,
                std::vector<GenomicRegion> &outputGRs) {
  outputGRs.clear();

  std::size_t covered_bases = 0;
  std::size_t read_iterator = inputMR.r.get_start();
  std::size_t seq_iterator = 0;
  std::size_t total_covered_bases = 0;

  while (seq_iterator < std::size(inputMR.seq)) {
    if (inputMR.seq[seq_iterator] != 'N')
      ++covered_bases;

    // if we reach the end of a bin, probabilistically create a binned read
    // with probability proportional to the number of covered bases
    if (read_iterator % bin_size == bin_size - 1) {
      const double frac = static_cast<double>(covered_bases) / bin_size;
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      if (dist(generator) <= frac) {
        const std::size_t curr_start =
          read_iterator - (read_iterator % bin_size);
        const std::size_t curr_end = curr_start + bin_size;
        const GenomicRegion binned_gr(
          inputMR.r.get_chrom(), curr_start, curr_end, inputMR.r.get_name(),
          inputMR.r.get_score(), inputMR.r.get_strand());
        outputGRs.push_back(binned_gr);
      }
      covered_bases = 0;
    }
    ++seq_iterator;
    ++read_iterator;
  }

  const double frac = static_cast<double>(covered_bases) / bin_size;
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  if (dist(generator) <= frac) {
    const std::size_t curr_start = read_iterator - (read_iterator % bin_size);
    const std::size_t curr_end = curr_start + bin_size;
    const GenomicRegion binned_gr(inputMR.r.get_chrom(), curr_start, curr_end,
                                  inputMR.r.get_name(), inputMR.r.get_score(),
                                  inputMR.r.get_strand());
    outputGRs.push_back(binned_gr);
  }
}

auto
load_coverage_counts_MR(const std::string &input_file_name,
                        const std::uint64_t seed, const std::size_t bin_size,
                        const std::size_t max_width,
                        std::vector<double> &coverage_hist) -> std::size_t {
  srand(time(nullptr) + getpid());
  // Runif runif(rand());
  std::mt19937 generator(seed);

  std::ifstream in(input_file_name);
  if (!in)
    throw std::runtime_error("problem opening file: " + input_file_name);

  MappedRead mr;
  if (!(in >> mr))
    throw std::runtime_error("problem reading from: " + input_file_name);

  // initialize prioirty queue to reorder the split reads
  ReadPQ PQ;

  std::size_t n_reads = 0;
  std::size_t n_bins = 0;
  GenomicRegion curr, prev;
  std::size_t current_count = 1;

  do {
    if (mr.r.get_width() > max_width)
      throw std::runtime_error("Encountered read of width " +
                               toa(mr.r.get_width()) +
                               "max_width set too small");

    std::vector<GenomicRegion> splitGRs;
    SplitMappedRead(mr, generator, bin_size, splitGRs);

    ++n_reads;
    n_bins += std::size(splitGRs);

    // add split Genomic Regions to the priority queue
    for (const auto &splitGR : splitGRs)
      PQ.push(splitGR);

    // remove Genomic Regions from the priority queue
    if (std::size(splitGRs) > 0)
      while (!PQ.empty() && is_ready_to_pop(PQ, splitGRs.back(), max_width))
        empty_pq(curr, prev, current_count, coverage_hist, PQ, input_file_name);
  } while (in >> mr);

  // done adding reads, now spit the rest out
  while (!PQ.empty())
    empty_pq(curr, prev, current_count, coverage_hist, PQ, input_file_name);

  return n_reads;
}

auto
load_coverage_counts_GR(const std::string &input_file_name,
                        const std::uint64_t seed, const std::size_t bin_size,
                        const std::size_t max_width,
                        std::vector<double> &coverage_hist) -> std::size_t {
  srand(time(nullptr) + getpid());
  // Runif runif(rand());
  std::mt19937 generator(seed);

  std::ifstream in(input_file_name);
  if (!in)
    throw std::runtime_error("problem opening file: " + input_file_name);

  GenomicRegion inputGR;
  if (!(in >> inputGR))
    throw std::runtime_error("problem reading from: " + input_file_name);

  // initialize prioirty queue to reorder the split reads
  ReadPQ PQ;

  // prev and current Genomic Regions to compare
  GenomicRegion curr;
  GenomicRegion prev;
  std::size_t n_reads = 0;
  std::size_t current_count = 1;

  do {
    std::vector<GenomicRegion> splitGRs;
    SplitGenomicRegion(inputGR, generator, bin_size, splitGRs);
    const auto n_splits = std::size(splitGRs);
    assert(!splitGRs.empty());
    const auto last_part = splitGRs.back();

    // add split Genomic Regions to the priority queue
    for (auto &&splitGR : splitGRs)
      PQ.push(std::move(splitGR));

    if (n_splits > 0) {  // ADS: is this a bug?
      // remove Genomic Regions from the priority queue
      while (!PQ.empty() && is_ready_to_pop(PQ, last_part, max_width))
        empty_pq(curr, prev, current_count, coverage_hist, PQ, input_file_name);
    }
    ++n_reads;
  } while (in >> inputGR);

  // done adding reads, now spit the rest out
  while (!PQ.empty())
    empty_pq(curr, prev, current_count, coverage_hist, PQ, input_file_name);

  return n_reads;
}

#ifdef HAVE_HTSLIB
// Deal with SAM/BAM format only if we have htslib

static inline bool
not_mapped(const bamxx::bam_rec &aln) {
  return get_tid(aln) == -1;
}

static inline void
swap(bamxx::bam_rec &a, bamxx::bam_rec &b) {
  std::swap(a.b, b.b);
}

struct aln_pos {
  int32_t tid{};
  hts_pos_t pos{};
  aln_pos() = default;
  aln_pos(const int32_t tid, const hts_pos_t pos) : tid{tid}, pos{pos} {}
  explicit aln_pos(const bamxx::bam_rec &a) :
    tid{get_tid(a)}, pos{get_pos(a)} {}
  bool operator<(const aln_pos &rhs) const {
    return tid < rhs.tid || (tid == rhs.tid && pos < rhs.pos);
  }
  bool operator>(const aln_pos &rhs) const {
    return tid > rhs.tid || (tid == rhs.tid && pos > rhs.pos);
  }
  bool operator!=(const aln_pos &rhs) const {
    // ADS: ordered to check pos first
    return pos != rhs.pos || tid != rhs.tid;
  }
};

struct aln_pos_pair {
  int32_t tid{};
  hts_pos_t pos{};
  int32_t mtid{};
  hts_pos_t mpos{};
  explicit aln_pos_pair(const bamxx::bam_rec &a) :
    tid{get_tid(a)}, pos{get_pos(a)}, mtid{get_mtid(a)}, mpos{get_mpos(a)} {}
  bool operator<(const aln_pos_pair &rhs) const {
    // ADS: only compares on tid and pos, NOT mtid or mpos
    return tid < rhs.tid || (tid == rhs.tid && pos < rhs.pos);
  }
  bool operator!=(const aln_pos_pair &rhs) const {
    // ADS: ordered to check pos first
    return pos != rhs.pos || tid != rhs.tid || mtid != rhs.mtid ||
           mpos != rhs.mpos;
  }
};

template <typename T>
static inline void
update_duplicate_counts_hist_BAM(const T &curr, const T &prev,
                                 vector<double> &counts_hist,
                                 size_t &current_count) {
  if (prev != curr) {
    // next read is new, update counts_hist to include current_count
    if (size(counts_hist) < current_count + 1) {
      // histogram is too small, resize
      counts_hist.resize(current_count + 1, 0.0);
    }
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

template <typename aln_pos_t>
size_t
load_counts_BAM(const uint32_t n_threads, const string &inputfile,
                vector<double> &counts_hist) {
  bamxx::bam_tpool tp(n_threads);

  bamxx::bam_in hts(inputfile);  // assume already checked
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw runtime_error("failed to read header");

  if (n_threads > 1)
    tp.set_io(hts);

  // find first mapped read to start
  bamxx::bam_rec aln;
  while (hts.read(hdr, aln) && not_mapped(aln))
    ;

  size_t n_reads{};
  // if all reads unmapped, must return
  if (not_mapped(aln))
    return n_reads;

  // to check that reads are sorted properly
  vector<bool> chroms_seen(get_n_targets(hdr), false);

  // start with prev_aln being first read
  aln_pos_t prev{aln};

  // start with count of 1 for first read seen
  size_t current_count = 1;

  while (hts.read(hdr, aln)) {
    if (not_mapped(aln))
      continue;  // skip unmapped reads

    const aln_pos_t curr{aln};

    // check that reads are sorted
    if (curr < prev)
      throw runtime_error("locations unsorted in: " + inputfile);

    if (curr.tid != prev.tid) {  // check that reads are sorted
      if (chroms_seen[curr.tid])
        throw runtime_error("input not sorted");
      chroms_seen[curr.tid] = true;
    }

    // check that mapped read is not secondary
    update_duplicate_counts_hist_BAM(curr, prev, counts_hist, current_count);
    ++n_reads;
    prev = curr;
  }

  // account for the last read
  if (size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return n_reads;
}

size_t
load_counts_BAM_se(const uint32_t n_threads, const string &inputfile,
                   vector<double> &counts_hist) {
  return load_counts_BAM<aln_pos>(n_threads, inputfile, counts_hist);
}

size_t
load_counts_BAM_pe(const uint32_t n_threads, const string &inputfile,
                   vector<double> &counts_hist) {
  return load_counts_BAM<aln_pos_pair>(n_threads, inputfile, counts_hist);
}

struct genomic_interval {
  int32_t tid{};  // indicates uninitialized
  hts_pos_t start{};
  hts_pos_t stop{};
  bool operator<(const genomic_interval &rhs) const {
    // clang-format off
    return (tid < rhs.tid ||
            (tid == rhs.tid &&
             (start < rhs.start ||
              (start == rhs.start &&
               (stop < rhs.stop)))));
    // clang-format on
  }
};

static inline uint32_t
size(const genomic_interval &gi) {
  return gi.stop - gi.start;
}

template <typename T>
static inline T
round_prob(const T x, const uint32_t bin_size, const double frac) {
  // probabilisticly round read ends so they are at bin boundaries
  const double lo = (x / bin_size) * bin_size;
  const double hi = ((x + bin_size - 1) / bin_size) * bin_size;
  return frac < (x - lo) ? lo : hi;
}

// split a mapped read into multiple genomic intervals based on the
// number of base pairs in each
static void
split_genomic_interval(const genomic_interval &gi, mt19937 &generator,
                       const hts_pos_t bin_size, vector<aln_pos> &output) {
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  // could shorten or lengthen; postcond: ends are at bin boundaries
  const hts_pos_t r_start = round_prob(gi.start, bin_size, dist(generator));
  const hts_pos_t r_stop = round_prob(gi.stop, bin_size, dist(generator));

  // gather all the parts at bin offsets
  for (auto pos = r_start; pos < r_stop; pos += bin_size)
    output.emplace_back(gi.tid, pos);
}

template <class T, class U>
static inline bool
can_pop(const T &pq, const U &last, const hts_pos_t max_dist) {
  return pq.top().tid != last.tid || pq.top().pos + max_dist < last.pos;
}

template <class T>
static void
update_coverage_hist(const T &curr, const T &prev, vector<double> &counts_hist,
                     size_t &current_count) {
  if (curr != prev) {
    if (counts_hist.size() < current_count + 1)  // histogram too small
      counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

// ADS: don't care if mapped reads are SE or PE, we only need the
// first mate for each mapped read
size_t
load_coverage_counts_BAM(const uint32_t n_threads, const string &inputfile,
                         const uint32_t seed, const size_t bin_size,
                         const size_t max_width,
                         vector<double> &coverage_hist) {
  std::mt19937 generator(seed);

  bamxx::bam_tpool tp(n_threads);
  bamxx::bam_in hts(inputfile);  // assume already checked
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw runtime_error("failed to read header");

  if (n_threads > 1)
    tp.set_io(hts);

  // find first mapped read to start
  bamxx::bam_rec aln;
  while (hts.read(hdr, aln) && not_mapped(aln))
    ;

  size_t n_reads{};
  if (not_mapped(aln))  // no reads unmapped
    return 0;

  // to check reads are sorted properly
  vector<bool> chroms_seen(get_n_targets(hdr), false);

  // start with count of 1 for first read seen
  size_t current_count = 1;

  // initialize prioirty queue to reorder the split reads
  priority_queue<aln_pos, vector<aln_pos>, std::greater<aln_pos>> pq;
  vector<aln_pos> parts;  // reuse allocated space
  aln_pos prev_part;
  genomic_interval prev;

  // max_dist indicates when we think we can assume the read parts
  // will be sorted and can be processed; this is not the same as the
  // full reads being sorted
  const hts_pos_t max_dist = bin_size + max_width;

  while (hts.read(hdr, aln)) {
    if (not_mapped(aln))
      continue;  // check that read is mapped

    const hts_pos_t len = rlen_from_cigar(aln);
    const genomic_interval curr{get_tid(aln), get_pos(aln), get_pos(aln) + len};

    if (curr.tid != prev.tid) {
      if (chroms_seen[curr.tid])
        throw runtime_error("input not sorted");
      chroms_seen[curr.tid] = true;
    }

    if (size(curr) > max_width)
      throw runtime_error("read " + string(bam_get_qname(aln)) + " covers " +
                          std::to_string(size(curr)) +
                          "bp; increase max width or reconsider data");

    parts.clear();  // need new vec, but keep capacity
    split_genomic_interval(curr, generator, bin_size, parts);

    // add split intervals to the priority queue
    const auto last = parts.back();  // keep a copy for test below
    for (const auto &i : parts)
      pq.push(i);

    // remove genomic interval parts from the priority queue
    while (!pq.empty() && can_pop(pq, last, max_dist)) {
      const aln_pos curr_part = pq.top();
      pq.pop();
      // update counts hist
      update_coverage_hist(curr_part, prev_part, coverage_hist, current_count);
      prev_part = curr_part;
    }
    prev = curr;
    ++n_reads;
  }

  // take care of remaining parts in priority queue
  while (!pq.empty()) {
    const aln_pos curr_part = pq.top();
    pq.pop();
    // update counts hist
    update_coverage_hist(curr_part, prev_part, coverage_hist, current_count);
    prev_part = curr_part;
  }
  return n_reads;
}

#endif  // HAVE_HTSLIB
