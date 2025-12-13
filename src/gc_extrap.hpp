/* Copyright (C) 2013-2025 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Authors: Timothy Daley and Andrew Smith
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_GC_EXTRAP_HPP_
#define SRC_GC_EXTRAP_HPP_

namespace gc_extrap {

static constexpr auto about_msg = R"(
Estimate the size of the part of the genome to be covered by mapped reads.
)";

static constexpr auto footer_msg = R"(
This approach is described in Daley & Smith (2014). The method is the same as
for lc_extrap: using rational function approximation to a power-series
expansion for the number of "unobserved" bases in the initial sample. The
gc_extrap method is adapted to deal with individual nucleotides rather than
distinct reads.
)";

auto
main(int argc, char *argv[]) -> int;

};  // namespace gc_extrap

#endif  // SRC_GC_EXTRAP_HPP_
