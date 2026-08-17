/* Copyright (C) 2013-2026 Andrew D. Smith
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

#ifndef SRC_LOAD_DATA_FOR_COMPLEXITY_HPP_
#define SRC_LOAD_DATA_FOR_COMPLEXITY_HPP_

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

enum class input_format_type : std::uint8_t {
  unknown,
  bam,
  bed,
  hist,
  counts,
};

[[nodiscard]] inline auto
is_unknown(input_format_type t) -> bool {
  return t == input_format_type::unknown;
}

[[nodiscard]] inline auto
is_bam(input_format_type t) -> bool {
  return t == input_format_type::bam;
}

[[nodiscard]] inline auto
is_bed(input_format_type t) -> bool {
  return t == input_format_type::bed;
}

[[nodiscard]] inline auto
is_hist(input_format_type t) -> bool {
  return t == input_format_type::hist;
}

[[nodiscard]] inline auto
is_counts(input_format_type t) -> bool {
  return t == input_format_type::counts;
}

[[nodiscard]] auto
to_string(input_format_type t) -> std::string;

[[nodiscard]] auto
get_input_format_type(const std::string &filename) -> input_format_type;

[[nodiscard]] auto
is_sam_or_bam_format(const std::string &filename) -> bool;

[[nodiscard]] auto
load_coverage_counts(const std::string &input_file_name,
                     const std::uint32_t seed, const std::size_t bin_size,
                     const std::size_t max_width)
  -> std::tuple<std::size_t, std::vector<double>>;

[[nodiscard]] auto
load_histogram(const std::string &filename)
  -> std::tuple<std::size_t, std::vector<double>>;

[[nodiscard]] auto
load_counts(const std::string &input_file_name)
  -> std::tuple<std::size_t, std::vector<double>>;

[[nodiscard]] auto
load_counts_bed_pe(const std::string &input_file_name)
  -> std::tuple<std::size_t, std::vector<double>>;

[[nodiscard]] auto
load_counts_bed_se(const std::string &input_file_name)
  -> std::tuple<std::size_t, std::vector<double>>;

#ifdef HAVE_HTSLIB

[[nodiscard]] auto
load_counts_BAM_pe(const std::uint32_t n_threads,
                   const std::string &input_file_name)
  -> std::tuple<std::size_t, std::vector<double>>;

[[nodiscard]] auto
load_counts_BAM_se(const std::uint32_t n_threads,
                   const std::string &input_file_name)
  -> std::tuple<std::size_t, std::vector<double>>;

[[nodiscard]] auto
load_coverage_counts_BAM(const std::uint32_t n_threads,
                         const std::string &input_file_name,
                         const std::uint32_t seed, const std::size_t bin_size,
                         const std::size_t max_width)
  -> std::tuple<std::size_t, std::vector<double>>;

#endif  // HAVE_HTSLIB

#endif  // SRC_LOAD_DATA_FOR_COMPLEXITY_HPP_
