/* preseq: to predict properties of genomic sequencing libraries
 *
 * Copyright (C) 2013-2024 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Authors: Timothy Daley, Chao Deng, Victoria Helus, and Andrew Smith
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

#include <config.h>

#include "common.hpp"

// the preseq commands
#include "bound_pop.hpp"
#include "c_curve.hpp"
#include "gc_extrap.hpp"
#include "lc_extrap.hpp"
#include "pop_size.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

static std::string
usage_message() {
  std::ostringstream oss;
  oss << "preseq: a program for analyzing library complexity\n"
         "Version: ";
  oss << VERSION;
  oss << "\n\n"
         "Usage: preseq <command> [OPTIONS]\n\n"
         "<command>: c_curve    generate complexity curve for a library\n"
         "           lc_extrap  predict the yield for future experiments\n"
         "           gc_extrap  predict genome coverage low input\n"
         "                      sequencing experiments\n"
         "           bound_pop  lower bound on population size\n"
         "           pop_size   estimate number of unique species\n";
  return oss.str();
}

int
main(const int argc, const char *argv[]) {
  if (argc < 2) {
    std::cerr << usage_message() << std::endl;
    return EXIT_SUCCESS;
  }

  static const std::string cmd = argv[1];

  if (cmd == "lc_extrap")
    return lc_extrap_main(argc, argv);

  if (cmd == "c_curve")
    return c_curve_main(argc, argv);

  if (cmd == "gc_extrap")
    return gc_extrap_main(argc, argv);

  if (cmd == "bound_pop")
    return bound_pop_main(argc, argv);

  if (cmd == "pop_size")
    return pop_size_main(argc, argv);

  std::cerr << "Error: unrecognized command: " << argv[1] << std::endl
            << usage_message() << std::endl;

  return EXIT_FAILURE;
}
