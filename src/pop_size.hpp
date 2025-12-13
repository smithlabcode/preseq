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

#ifndef SRC_POP_SIZE_HPP_
#define SRC_POP_SIZE_HPP_

namespace pop_size {

static constexpr auto about_msg = R"(
Estimate the population size using a small sample from the population.
)";

static constexpr auto footer_msg = R"(
Estimate the total population size using the approach described in Daley &
Smith (2013), extrapolating to very long range. Default parameters assume that
the initial sample represents at least 1e-9 of the population, which is
sufficient for every example application we have seen.
)";

auto
main(int argc, char *argv[]) -> int;

};  // namespace pop_size

#endif  // SRC_POP_SIZE_HPP_
