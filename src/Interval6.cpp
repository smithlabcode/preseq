/* Copyright (C) 2025 Andrew D Smith
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this software; if not, write to the Free Software Foundation, Inc., 51
 * Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "Interval6.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

auto
Interval6::initialize(const char *c, const char *c_end) -> bool {
  constexpr auto is_sep = [](const char x) { return x == ' ' || x == '\t'; };
  constexpr auto not_sep = [](const char x) { return x != ' ' && x != '\t'; };

  bool failed = false;

  // NOLINTBEGIN(*-pointer-arithmetic)
  auto field_s = c;
  auto field_e = std::find_if(field_s + 1, c_end, is_sep);
  if (field_e == c_end)
    failed = true;

  // chrom
  {
    const std::uint32_t d = std::distance(field_s, field_e);
    chrom = std::string{field_s, d};
  }

  // start
  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);
  {
    const auto [ptr, ec] = std::from_chars(field_s, field_e, start);
    failed = failed || (ptr == field_s);
  }

  // stop
  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);
  {
    const auto [ptr, ec] = std::from_chars(field_s, field_e, stop);
    failed = failed || (ptr == field_s);
  }

  // name
  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);
  name = std::string(field_s, std::distance(field_s, field_e));

  // score
  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);
  {
#ifdef __APPLE__
    const int ret = std::sscanf(field_s, "%lf", &score);
    failed = failed || (ret < 1);
#else
    const auto [ptr, ec] = std::from_chars(field_s, field_e, score);
    failed = failed || ec != std::errc{};
#endif
  }

  // strand (no stop; just one char and maybe end of line)
  field_s = std::find_if(field_e + 1, c_end, not_sep);
  failed = failed || (field_s == c_end);
  strand = *field_s;
  failed = failed || (strand != '-' && strand != '+');
  // NOLINTEND(*-pointer-arithmetic)

  return !failed;
}

[[nodiscard]] auto
read_intervals6(const std::string &intervals_file) -> std::vector<Interval6> {
  std::ifstream in(intervals_file);
  if (!in)
    throw std::runtime_error("failed to open file: " + intervals_file);
  std::string line;
  std::vector<Interval6> intervals;
  while (getline(in, line))
    intervals.emplace_back(line);
  return intervals;
}
