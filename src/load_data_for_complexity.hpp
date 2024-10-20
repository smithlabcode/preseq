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

#ifndef SRC_LOAD_DATA_FOR_COMPLEXITY_HPP_
#define SRC_LOAD_DATA_FOR_COMPLEXITY_HPP_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

std::size_t
load_coverage_counts_MR(const std::string &input_file_name,
                        const std::uint32_t seed, const std::size_t bin_size,
                        const std::size_t max_width,
                        std::vector<double> &coverage_hist);

std::size_t
load_coverage_counts_GR(const std::string &input_file_name,
                        const std::uint32_t seed, const std::size_t bin_size,
                        const std::size_t max_width,
                        std::vector<double> &coverage_hist);

std::size_t
load_histogram(const std::string &filename, std::vector<double> &counts_hist);

std::size_t
load_counts(const std::string &input_file_name,
            std::vector<double> &counts_hist);

std::size_t
load_counts_BED_pe(const std::string &input_file_name,
                   std::vector<double> &counts_hist);

std::size_t
load_counts_BED_se(const std::string &input_file_name,
                   std::vector<double> &counts_hist);

#ifdef HAVE_HTSLIB
std::size_t
load_counts_BAM_pe(const std::uint32_t n_threads,
                   const std::string &input_file_name,
                   std::vector<double> &counts_hist);

std::size_t
load_counts_BAM_se(const std::uint32_t n_threads,
                   const std::string &input_file_name,
                   std::vector<double> &counts_hist);

std::size_t
load_coverage_counts_BAM(const std::uint32_t n_threads,
                         const std::string &input_file_name,
                         const std::uint32_t seed, const std::size_t bin_size,
                         const std::size_t max_width,
                         std::vector<double> &coverage_hist);

#endif  // HAVE_HTSLIB

#endif  // SRC_LOAD_DATA_FOR_COMPLEXITY_HPP_
