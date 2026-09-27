// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "bound_pop.hpp"
#include "c_curve.hpp"
#include "cli_common.hpp"
#include "common.hpp"
#include "gc_extrap.hpp"
#include "lc_extrap.hpp"
#include "pop_size.hpp"

#include <config.h>

#include <CLI11/CLI11.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <print>
#include <string>

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
  app.footer(description);
  // if (argc >= 3)
  //   app.footer(footer_msg);

  app.require_subcommand(0, 1);
  app.allow_extras();

  bool print_version{};

  // clang-format off
  app.add_flag("--version", print_version, "output version information and exit");
#ifdef INCLUDE_FULL_LICENSE_INFO
  app.add_flag("--licenses", print_licenses, "view licenses");
#endif
  const auto lc_extrap_cmd = app.add_subcommand("lc_extrap", lc_extrap_about_msg);
  const auto gc_extrap_cmd = app.add_subcommand("gc_extrap", gc_extrap_about_msg);
  const auto pop_size_cmd = app.add_subcommand("pop_size", pop_size_about_msg);
  const auto bound_pop_cmd = app.add_subcommand("bound_pop", bound_pop_about_msg);
  const auto c_curve_cmd = app.add_subcommand("c_curve", c_curve_about_msg);
  // clang-format on

  if (argc < 2) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  if (print_version) {
    std::println("{}", VERSION);
    return EXIT_SUCCESS;
  }

#ifdef INCLUDE_FULL_LICENSE_INFO
  if (view_licenses) {
    std::println("{}", license_text);
    return EXIT_SUCCESS;
  }
#endif

  if (app.got_subcommand(lc_extrap_cmd))
    return lc_extrap_main(argc - 1, argv + 1);

  if (app.got_subcommand(gc_extrap_cmd))
    return gc_extrap_main(argc - 1, argv + 1);

  if (app.got_subcommand(c_curve_cmd))
    return c_curve_main(argc - 1, argv + 1);

  if (app.got_subcommand(pop_size_cmd))
    return pop_size_main(argc - 1, argv + 1);

  if (app.got_subcommand(bound_pop_cmd))
    return bound_pop_main(argc - 1, argv + 1);

  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  std::println(std::cerr, "unrecognized command: {}", argv[1]);
  return EXIT_SUCCESS;
}
