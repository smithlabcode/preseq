/* preseq: a tool for analyzing sequencing library complexity
 *
 * Copyright (C) 2013-2025 University of Southern California and
 *                         Andrew D. Smith and Timothy Daley
 *
 * Author: Andrew D Smith
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

#include "bound_pop.hpp"
#include "c_curve.hpp"
#include "common.hpp"
#include "gc_extrap.hpp"
#include "lc_extrap.hpp"
#include "pop_size.hpp"

#include "CLI11/CLI11.hpp"

#include <config.h>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#ifdef INCLUDE_FULL_LICENSE_INFO
#include <license.h>
#endif

const auto description = R"(
Extrapolate the complexity of a library. This is the approach described in
Daley & Smith (2013). The method applies rational function approximation via
continued fractions with the original goal of estimating the number of
distinct reads that a sequencing library would yield upon deeper sequencing.
This method has been used for many different purposes since then.
)";

int
main(int argc, char *argv[]) {  // NOLINT(*-c-arrays)
  CLI::App app{"preseq: a tool for analyzing sequencing library complexity"};
  argv = app.ensure_utf8(argv);
  app.formatter(std::make_shared<preseq_formatter>());
  app.usage("\nUsage: preseq command [OPTIONS]");
  // if (argc >= 2)
  // app.footer(rlstrip(description));
  // if (argc >= 3)
  //   app.footer(rlstrip(footer_msg));

  app.require_subcommand(0, 1);
  app.allow_extras();

  bool print_version{};

  // clang-format off
  app.add_flag("--version", print_version, "output version information and exit");
#ifdef INCLUDE_FULL_LICENSE_INFO
  app.add_flag("--licenses", print_licenses, "view licenses");
#endif
  const auto lc_extrap = app.add_subcommand("lc_extrap", rlstrip(lc_extrap::about_msg));
  const auto gc_extrap = app.add_subcommand("gc_extrap", rlstrip(gc_extrap::about_msg));
  const auto pop_size = app.add_subcommand("pop_size", rlstrip(pop_size::about_msg));
  const auto bound_pop = app.add_subcommand("bound_pop", rlstrip(bound_pop::about_msg));
  const auto c_curve = app.add_subcommand("c_curve", rlstrip(c_curve::about_msg));
  // clang-format on

  if (argc < 2) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  if (print_version) {
    std::cout << VERSION << '\n';
    return EXIT_SUCCESS;
  }

#ifdef INCLUDE_FULL_LICENSE_INFO
  if (view_licenses) {
    std::cout << license_text;
    return EXIT_SUCCESS;
  }
#endif

  if (app.got_subcommand(lc_extrap))
    return lc_extrap::main(argc - 1, argv + 1);

  if (app.got_subcommand(gc_extrap))
    return gc_extrap::main(argc - 1, argv + 1);

  if (app.got_subcommand(c_curve))
    return c_curve::main(argc - 1, argv + 1);

  if (app.got_subcommand(pop_size))
    return pop_size::main(argc - 1, argv + 1);

  std::println(std::cerr, "unrecognized command: {}", std::string(argv[1]));
  std::println(std::cerr, "{}", usage_message);
  return EXIT_SUCCESS;
}
