/* Copyright (C) 2014-2025 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Authors: Andrew D. Smith and Timothy Daley
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "load_data_for_complexity.hpp"

#include "Interval6.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>  // IWYU pragma: keep
#include <iterator>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>  // IWYU pragma: keep
#include <utility>
#include <vector>

#ifdef HAVE_HTSLIB
#include "bam_record_utils.hpp"
#include <bamxx.hpp>
#include <htslib/sam.h>
#endif

// NOLINTBEGIN(*-narrowing-conversions)

[[nodiscard]] auto
to_string(input_format_type t) -> std::string {
  switch (t) {
  case input_format_type::bam:
    return "BAM";
  case input_format_type::bed:
    return "BED";
  case input_format_type::hist:
    return "HIST";
  case input_format_type::counts:
    return "COUNTS";
  default:  // case input_format_type::counts:
    return "UNKNOWN";
  }
  return std::string();
}

[[nodiscard]] auto
is_sam_or_bam_format(const std::string &filename) -> bool {
  // make sure the file can be opened at all
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("failed to open file: " + filename);
#ifdef HAVE_HTSLIB
  // if the file can be opened and has a header, we can get reads from it
  bamxx::bam_in hts(filename);
  if (!hts)
    return false;
  bamxx::bam_header hdr(hts);
  if (hdr)
    return true;
#endif
  return false;
}

[[nodiscard]] auto
get_input_format_type(const std::string &filename) -> input_format_type {
  if (is_sam_or_bam_format(filename))
    return input_format_type::bam;

  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("failed to open file: " + filename);

  std::string line;
  if (!std::getline(in, line))
    return input_format_type::unknown;
  std::string s;
  std::uint32_t val1{}, val2{};
  {
    std::istringstream iss(line);
    if (iss >> s >> val1 >> val2)
      return input_format_type::bed;
  }
  {
    std::istringstream iss(line);
    if (iss >> val1 >> val2)
      return input_format_type::hist;
  }
  {
    std::istringstream iss(line);
    if (iss >> val1)
      return input_format_type::bed;
  }
  return input_format_type::unknown;
}

template <typename T>
[[nodiscard]] static inline auto
width(const T &x) -> std::uint32_t {
  return x.stop - x.start;
}

static auto
update_pe_duplicate_counts_hist(const Interval6 &curr, const Interval6 &prev,
                                std::vector<double> &counts_hist,
                                std::size_t &current_count) -> bool {
  // check if reads are sorted
  if (curr.chrom == prev.chrom && curr.start < prev.start &&
      curr.stop < prev.stop)
    return false;

  // check if next read is new, and if so update counts_hist to include
  // current_count
  if (curr.chrom != prev.chrom || curr.start != prev.start ||
      curr.stop != prev.stop) {
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
update_se_duplicate_counts_hist(const Interval6 &curr, const Interval6 &prev,
                                const std::string &input_file_name,
                                std::vector<double> &counts_hist,
                                std::size_t &current_count) {
  // check if reads are sorted
  if (curr.chrom == prev.chrom && curr.start < prev.start)
    throw std::runtime_error("locations unsorted in: " + input_file_name);

  if (curr.chrom != prev.chrom || curr.start != prev.start) {
    // next read is new, update counts_hist to include current_count
    // histogram is too small, resize
    if (std::size(counts_hist) < current_count + 1)
      counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

struct interval_greater {
  auto
  operator()(const Interval6 &a, const Interval6 &b) const -> bool {
    return b < a;  // i.e. a > b
  }
};

using read_pq =
  std::priority_queue<Interval6, std::vector<Interval6>, interval_greater>;

static auto
is_ready_to_pop(const read_pq &pq, const Interval6 &interval,
                const std::size_t max_width) -> bool {
  return pq.top().chrom != interval.chrom ||
         pq.top().stop + max_width < interval.start;
}

static void
empty_pq(Interval6 &prev, std::size_t &current_count,
         std::vector<double> &counts_hist, read_pq &read_pq,
         const std::string &input_file_name) {
  const auto curr = read_pq.top();
  read_pq.pop();

  // update counts hist
  const bool update_success =
    update_pe_duplicate_counts_hist(curr, prev, counts_hist, current_count);
  if (!update_success) {
    std::ostringstream oss;
    oss << "reads unsorted in: " << input_file_name << "\n"
        << "prev = \t" << to_string(prev) << "\n"
        << "curr = \t" << to_string(curr) << "\n"
        << "Increase seg_len if in paired end mode";
    throw std::runtime_error(oss.str());
  }
  prev = curr;
}

// for BED file input

[[nodiscard]] auto
load_counts_bed_se(const std::string &input_file_name)
  -> std::tuple<std::size_t, std::vector<double>> {
  std::vector<double> counts_hist(2, 0.0);

  std::ifstream in(input_file_name);
  if (!in)
    throw std::runtime_error("problem opening file: " + input_file_name);

  std::size_t n_reads{};
  std::size_t current_count{};

  Interval6 prev;
  std::string line;
  while (std::getline(in, line)) {
    const auto curr = Interval6(line);
    update_se_duplicate_counts_hist(curr, prev, input_file_name, counts_hist,
                                    current_count);
    ++n_reads;
    prev = curr;
  }

  // to account for the last read compared to the one before it.
  if (std::size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return std::make_tuple(n_reads, std::move(counts_hist));
}

auto
load_counts_bed_pe(const std::string &input_file_name)
  -> std::tuple<std::size_t, std::vector<double>> {
  std::vector<double> counts_hist(2, 0.0);

  std::ifstream in(input_file_name);
  if (!in)
    throw std::runtime_error("problem opening file: " + input_file_name);

  std::size_t n_reads{};
  std::size_t current_count{};

  Interval6 prev;
  std::string line;

  // read in file and compare each gr with the one before it
  while (std::getline(in, line)) {
    const auto curr = Interval6(line);
    const bool update_success =
      update_pe_duplicate_counts_hist(curr, prev, counts_hist, current_count);
    if (!update_success)
      throw std::runtime_error("reads unsorted in " + input_file_name);
    ++n_reads;
    prev = curr;
  }

  if (std::size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);

  // to account for the last read compared to the one before it.
  ++counts_hist[current_count];

  return std::make_tuple(n_reads, std::move(counts_hist));
}

[[nodiscard]] auto
load_counts(const std::string &infile)
  -> std::tuple<std::size_t, std::vector<double>> {
  std::ifstream in(infile);
  if (!in)
    throw std::runtime_error("failed to open file: " + infile);

  std::vector<double> vals((std::istream_iterator<double>(in)),
                           std::istream_iterator<double>());
  if (vals.empty())
    return std::make_tuple(0, std::vector<double>{});

  const auto max_val = *std::max_element(std::cbegin(vals), std::cend(vals));

  std::vector<double> counts_hist(max_val + 1, 0.0);
  for (const auto v : vals)
    ++counts_hist[v];

  const auto n_reads = std::accumulate(std::cbegin(vals), std::cend(vals), 0ul);
  return std::make_tuple(n_reads, std::move(counts_hist));
}

// returns number of reads from file containing counts histogram
[[nodiscard]] auto
load_histogram(const std::string &filename)
  -> std::tuple<std::size_t, std::vector<double>> {
  std::ifstream in(filename);
  if (!in)  // if file doesn't open
    throw std::runtime_error("could not open histogram: " + filename);

  std::vector<double> counts_hist;

  std::size_t n_reads = 0;
  std::size_t line_count = 0ul, prev_read_count = 0ul;
  std::string buffer;
  while (std::getline(in, buffer)) {
    ++line_count;
    std::size_t read_count = 0ul;
    double frequency = 0.0;
    std::istringstream is(buffer);
    // error reading input
    if (!(is >> read_count >> frequency))
      throw std::runtime_error("bad histogram line format:\n" + buffer + "\n" +
                               "(line " + std::to_string(line_count) + ")");

    // histogram is out of order?
    if (read_count < prev_read_count)
      throw std::runtime_error("bad line order in file " + filename + "\n" +
                               "(line " + std::to_string(line_count) + ")");
    counts_hist.resize(read_count + 1, 0.0);
    counts_hist[read_count] = frequency;
    if (read_count == 0ul) {
      throw std::runtime_error("counts histograms may not "
                               "include an entry for zero");
    }
    prev_read_count = read_count;
    n_reads += static_cast<std::size_t>(read_count * frequency);
  }

  return std::make_tuple(n_reads, std::move(counts_hist));
}

// Loading coverage counts

// probabilistically split intervals into mutiple intervals of width
// equal to bin_size
[[nodiscard]] static auto
split_genomic_region(Interval6 interval, std::mt19937 &generator,
                     const std::uint32_t bin_size) -> std::vector<Interval6> {
  const auto frac = static_cast<double>(interval.start % bin_size) / bin_size;
  const auto w = width(interval);

  // ADS: this seems like a bunch of duplicated code just for a single
  // function difference
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  const auto start_diff = dist(generator) > frac ? 0 : bin_size - 1;
  interval.start = ((interval.start + start_diff) / bin_size) * bin_size;
  interval.stop = interval.start + w;

  std::vector<Interval6> parts;
  for (auto i = 0u; i < width(interval); i += bin_size) {
    const std::uint32_t curr_start = interval.start + i;
    const double curr_end = std::min(interval.stop, curr_start + bin_size);
    if (dist(generator) <= (curr_end - curr_start) / bin_size)
      parts.emplace_back(interval.chrom, curr_start, curr_start + bin_size,
                         interval.name, interval.score, interval.strand);
  }
  return parts;
}

[[nodiscard]] auto
load_coverage_counts(const std::string &infile, const std::uint32_t seed,
                     const std::size_t bin_size, const std::size_t max_width)
  -> std::tuple<std::size_t, std::vector<double>> {
  std::mt19937 generator(seed);

  std::ifstream in(infile);
  if (!in)
    throw std::runtime_error("problem opening file: " + infile);

  // prioirty queue to reorder the split reads
  read_pq pq;

  Interval6 prev;
  std::size_t n_reads{};
  std::size_t current_count{1};

  std::vector<double> coverage_hist;

  std::string line;
  while (std::getline(in, line)) {
    const auto interval = Interval6(line);
    const auto splits = split_genomic_region(interval, generator, bin_size);

    // add split intervals to the priority queue
    for (const auto &i : splits)
      pq.push(i);

    if (std::size(splits) > 0) {
      // remove intervals from the priority queue
      while (!pq.empty() && is_ready_to_pop(pq, splits.back(), max_width))
        empty_pq(prev, current_count, coverage_hist, pq, infile);
    }
    ++n_reads;
  }

  // done adding reads, now spit the rest out
  while (!pq.empty())
    empty_pq(prev, current_count, coverage_hist, pq, infile);

  return std::make_tuple(n_reads, std::move(coverage_hist));
}

#ifdef HAVE_HTSLIB

struct genomic_interval {
  std::int32_t tid{-1};  // indicates uninitialized
  hts_pos_t start{};
  hts_pos_t stop{};
  auto
  operator<(const genomic_interval &rhs) const -> bool {
    // clang-format off
    return (tid < rhs.tid ||
            (tid == rhs.tid &&
             (start < rhs.start ||
              (start == rhs.start &&
               (stop < rhs.stop)))));
    // clang-format on
  }
};

struct aln_pos {
  std::int32_t tid{};
  hts_pos_t pos{};
  aln_pos() = default;
  aln_pos(const std::int32_t tid, const hts_pos_t pos) : tid{tid}, pos{pos} {}
  explicit aln_pos(const bamxx::bam_rec &a) :
    tid{get_tid(a)}, pos{get_pos(a)} {}
  auto
  operator<(const aln_pos &rhs) const -> bool {
    return tid < rhs.tid || (tid == rhs.tid && pos < rhs.pos);
  }
  auto
  operator>(const aln_pos &rhs) const -> bool {
    return tid > rhs.tid || (tid == rhs.tid && pos > rhs.pos);
  }
  auto
  operator!=(const aln_pos &rhs) const -> bool {
    // ADS: ordered to check pos first
    return pos != rhs.pos || tid != rhs.tid;
  }
};

struct aln_pos_pair {
  std::int32_t tid{};
  hts_pos_t pos{};
  std::int32_t mtid{};
  hts_pos_t mpos{};
  explicit aln_pos_pair(const bamxx::bam_rec &a) :
    tid{get_tid(a)}, pos{get_pos(a)}, mtid{get_mtid(a)}, mpos{get_mpos(a)} {}
  auto
  operator<(const aln_pos_pair &rhs) const -> bool {
    // ADS: only compares on tid and pos, NOT mtid or mpos
    return tid < rhs.tid || (tid == rhs.tid && pos < rhs.pos);
  }
  auto
  operator!=(const aln_pos_pair &rhs) const -> bool {
    // ADS: ordered to check pos first
    return pos != rhs.pos || tid != rhs.tid || mtid != rhs.mtid ||
           mpos != rhs.mpos;
  }
};

template <typename T>
[[nodiscard]] static inline auto
round_position(const T x, const std::uint32_t bin_size,
               const double frac) -> T {
  // probabilisticly round read ends so they are at bin boundaries
  const double lo = (x / bin_size) * bin_size;
  const double hi = ((x + bin_size - 1) / bin_size) * bin_size;
  return frac < (x - lo) ? lo : hi;
}

// split a mapped read into multiple genomic intervals based on the number of
// base pairs in each
static auto
split_genomic_interval(const genomic_interval &gi, std::mt19937 &generator,
                       const hts_pos_t bin_size) -> std::vector<aln_pos> {
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  // could shorten or lengthen; postcond: ends are at bin boundaries
  const hts_pos_t r_start = round_position(gi.start, bin_size, dist(generator));
  const hts_pos_t r_stop = round_position(gi.stop, bin_size, dist(generator));

  // gather all the parts at bin offsets
  std::vector<aln_pos> parts;
  for (auto pos = r_start; pos < r_stop; pos += bin_size)
    parts.emplace_back(gi.tid, pos);
  return parts;
}

static inline auto
not_mapped(const bamxx::bam_rec &aln) -> bool {
  return get_tid(aln) == -1;
}

template <typename T>
static inline void
update_duplicate_counts_hist_BAM(const T &curr, const T &prev,
                                 std::vector<double> &counts_hist,
                                 std::size_t &current_count) {
  if (prev != curr) {
    // next read is new, update counts_hist to include current_count
    if (std::size(counts_hist) < current_count + 1) {
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
auto
load_counts_BAM(const std::uint32_t n_threads, const std::string &inputfile)
  -> std::tuple<std::size_t, std::vector<double>> {
  bamxx::bam_tpool tp(n_threads);

  bamxx::bam_in hts(inputfile);  // assume already checked
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  if (n_threads > 1)
    tp.set_io(hts);

  // find first mapped read to start
  bamxx::bam_rec aln;
  while (hts.read(hdr, aln) && not_mapped(aln))
    ;

  // if all reads unmapped, must return
  if (not_mapped(aln))
    return std::make_tuple(0, std::vector<double>{});

  // to check that reads are sorted properly
  std::vector<bool> chroms_seen(get_n_targets(hdr), false);

  // start with prev_aln being first read
  aln_pos_t prev{aln};

  // start with count of 1 for first read seen
  std::size_t n_reads{};
  std::size_t current_count{1};

  std::vector<double> counts_hist;

  while (hts.read(hdr, aln)) {
    if (not_mapped(aln))
      continue;  // skip unmapped reads

    const aln_pos_t curr{aln};

    // check that reads are sorted
    if (curr < prev)
      throw std::runtime_error("locations unsorted in: " + inputfile);

    if (curr.tid != prev.tid) {  // check that reads are sorted
      if (chroms_seen[curr.tid])
        throw std::runtime_error("input not sorted");
      chroms_seen[curr.tid] = true;
    }

    // check that mapped read is not secondary
    update_duplicate_counts_hist_BAM(curr, prev, counts_hist, current_count);
    ++n_reads;
    prev = curr;
  }

  // account for the last read
  if (std::size(counts_hist) < current_count + 1)
    counts_hist.resize(current_count + 1, 0.0);
  ++counts_hist[current_count];

  return std::make_tuple(n_reads, std::move(counts_hist));
}

[[nodiscard]] auto
load_counts_BAM_se(const std::uint32_t n_threads, const std::string &inputfile)
  -> std::tuple<std::size_t, std::vector<double>> {
  return load_counts_BAM<aln_pos>(n_threads, inputfile);
}

[[nodiscard]] auto
load_counts_BAM_pe(const std::uint32_t n_threads, const std::string &inputfile)
  -> std::tuple<std::size_t, std::vector<double>> {
  return load_counts_BAM<aln_pos_pair>(n_threads, inputfile);
}

template <class T>
static void
update_coverage_hist(const T &curr, const T &prev,
                     std::vector<double> &counts_hist,
                     std::size_t &current_count) {
  if (curr != prev) {
    if (std::size(counts_hist) < current_count + 1)  // histogram too small
      counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    current_count = 1;
  }
  else  // next read is same, update current_count
    ++current_count;
}

// ADS: don't care if mapped reads are SE or PE, we only need the first mate
// for each mapped read
auto
load_coverage_counts_BAM(
  const std::uint32_t n_threads, const std::string &inputfile,
  const std::uint32_t seed, const std::size_t bin_size,
  const std::size_t max_width) -> std::tuple<std::size_t, std::vector<double>> {
  std::mt19937 generator(seed);

  bamxx::bam_tpool tp(n_threads);
  bamxx::bam_in hts(inputfile);
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  if (n_threads > 1)
    tp.set_io(hts);

  // find first mapped read to start
  bamxx::bam_rec aln;
  while (hts.read(hdr, aln) && not_mapped(aln))
    ;

  std::size_t n_reads{};
  if (not_mapped(aln))  // no reads unmapped
    return {0, {}};

  // to check reads are sorted properly
  std::vector<bool> chroms_seen(get_n_targets(hdr), false);

  // start with count of 1 for first read seen
  std::size_t current_count = 1;

  // initialize prioirty queue to reorder the split reads
  std::priority_queue<aln_pos, std::vector<aln_pos>, std::greater<>> pq;
  aln_pos prev_part;
  genomic_interval prev;

  // max_dist indicates when we think we can assume the read parts will be
  // sorted and can be processed; not the same as the full reads being sorted
  const hts_pos_t max_dist = bin_size + max_width;

  std::vector<double> coverage_hist;

  while (hts.read(hdr, aln)) {
    if (not_mapped(aln))
      continue;  // check that read is mapped

    const hts_pos_t len = rlen_from_cigar(aln);
    const genomic_interval curr{get_tid(aln), get_pos(aln), get_pos(aln) + len};

    if (curr.tid != prev.tid) {
      if (chroms_seen[curr.tid])
        throw std::runtime_error("input not sorted");
      chroms_seen[curr.tid] = true;
    }

    if (width(curr) > max_width)
      throw std::runtime_error("read " + std::string(bam_get_qname(aln)) +
                               " covers " + std::to_string(width(curr)) +
                               "bp; increase max width or reconsider data");

    const auto parts = split_genomic_interval(curr, generator, bin_size);

    // add split intervals to the priority queue
    const auto last = parts.back();  // keep a copy for test below
    for (const auto &i : parts)
      pq.push(i);

    const auto can_pop = [&pq, &last, max_dist]() {
      return pq.top().tid != last.tid || pq.top().pos + max_dist < last.pos;
    };

    // remove genomic interval parts from the priority queue
    while (!pq.empty() && can_pop()) {
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
  return std::make_tuple(n_reads, std::move(coverage_hist));
}

#endif  // HAVE_HTSLIB

// NOLINTEND(*-narrowing-conversions)
