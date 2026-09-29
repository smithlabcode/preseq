// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "bound_pop.hpp"
#include "c_curve.hpp"
#include "cli_common.hpp"
#include "common.hpp"
#include "gc_extrap.hpp"
#include "lc_extrap.hpp"
#include "pop_size.hpp"

#include <config.h>
#include <license.h>

#include <CLI11/CLI11.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <print>
#include <span>
#include <string>

const auto description =
  R"(Estimate complexity characteristics of a DNA sequencing library.
Preseq includes several commands. The lc_extrap command is the most
general and can be used for applications outside of DNA sequencing.
)";

int
main(int argc, char *argv[]) {  // NOLINT(*-c-arrays)
  try {

    const std::span args(argv, argc);

    CLI::App app{"preseq: analyze DNA sequencing library complexity"};
    argv = app.ensure_utf8(argv);
    app.usage("Usage: preseq command [OPTIONS]");
    if (argc >= 2)
      app.footer(description);

    // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
    app.formatter(std::make_shared<preseq_formatter>());
    app.get_formatter()->column_width(16);
    app.get_formatter()->long_option_alignment_ratio(0.2);
    // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

    app.require_subcommand(0, 1);
    app.allow_extras();
    app.set_help_flag("");
    // clang-format off
    app.add_subcommand("lc_extrap", lc_extrap_about_msg)
      ->callback([&]{lc_extrap_main(args.subspan(1));});
    app.add_subcommand("c_curve", c_curve_about_msg)
      ->callback([&]{c_curve_main(args.subspan(1));});
    app.add_subcommand("gc_extrap", gc_extrap_about_msg)
      ->callback([&]{gc_extrap_main(args.subspan(1));});
    app.add_subcommand("pop_size", pop_size_about_msg)
      ->callback([&]{pop_size_main(args.subspan(1));});
    app.add_subcommand("bound_pop", bound_pop_about_msg)
      ->callback([&]{bound_pop_main(args.subspan(1));});
    // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
    app.set_version_flag("--version", VERSION, "Print program version");
    app.add_flag("--license", [&](auto) {
    std::print("{}", license_text); throw CLI::Success(); },
      "Print full license")
      ->callback_priority(CLI::CallbackPriority::PreRequirementsCheck);
    // clang-format on

    CLI11_PARSE(app, std::ssize(args), std::data(args));
    if (std::ssize(args) == 1) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }

    if (!app.get_subcommands().empty())
      return EXIT_SUCCESS;

    std::println("unrecognized command: {}", args[1]);
    return EXIT_FAILURE;
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
