/* Copyright (C) 2013-2024 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Authors: Timothy Daley and Andrew Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_LC_EXTRAP_HPP_
#define SRC_LC_EXTRAP_HPP_

constexpr auto lc_extrap_about_msg = R"(
Estimate a complexity curve (e.g., for a sequencing library).
)";

constexpr auto lc_extrap_footer_msg = R"(
This is the approach described in Daley & Smith (2013). The method applies
rational function approximation via continued fractions with the original goal
of estimating the number of distinct reads that a sequencing library would
yield upon deeper sequencing.  This method has been used for many different
purposes since then.
)";

auto
lc_extrap_main(int argc, char *argv[]) -> int;

#endif  // SRC_LC_EXTRAP_HPP_
