// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#include "c_curve.hpp"

#include "cli_common.hpp"
#include "common.hpp"
#include "load_data_for_complexity.hpp"

#include <CLI11/CLI11.hpp>
#include <libpreseq.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

auto
c_curve_main(const std::span<char *> args) -> int {
  static constexpr auto cmd_name = "c_curve";
  static constexpr auto n_points_default = 100;

  const int argc = std::ssize(args);
  auto argv = std::data(args);

  double step_size{};
  std::string outfile;
  std::string infile;
  bool verbose{};

  std::uint32_t n_threads{1};
  std::uint32_t n_points{n_points_default};

  CLI::App app{c_curve_about_msg};
  argv = app.ensure_utf8(argv);
  app.failure_message(
    [&](const CLI::App *a, const CLI::Error &b) -> std::string {
      return std::format("preseq command: {}\n", cmd_name) +
             CLI::FailureMessage::simple(a, b);
    });
  app.usage(std::format("Usage: preseq {} [OPTIONS]", cmd_name));
  if (argc >= 2)
    app.footer(" ");  // DESCRIPTION

  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  app.formatter(std::make_shared<preseq_formatter>());
  app.get_formatter()->column_width(22);
  app.get_formatter()->long_option_alignment_ratio(0.2);
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  // clang-format off
  app.add_option("INPUT", infile, "input file (hist/BED/BAM/SAM/values)")
    ->required()
    ->option_text("FILE")
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output", outfile, "Output file. Two columns, has header.")
    ->required()
    ->option_text("FILE");
  const auto step_size_opt =
    app.add_option("-s,--step", step_size,
                   "Step size. Must be smaller than data size.")
    ->option_text("FLOAT");
  app.add_option("-p,--points", n_points,
                 "Number of interpolation points. Excludes step size.")
    ->option_text("INT")
    ->excludes(step_size_opt);
  app.add_flag("-v,--verbose", verbose, "print more info about the data");
  // clang-format on

  if (argc == 1) {
    std::println("{}", app.help());
    return EXIT_SUCCESS;
  }
  CLI11_PARSE(app, argc, argv);

  const auto input_format = get_input_format_type(infile);
  if (is_unknown(input_format)) {
    std::println("unknown input format");
    return EXIT_FAILURE;
  }

  const auto [n_reads, counts_hist] = [&] {
    switch (input_format) {
    case input_format_type::hist:
      return load_histogram(infile);
    case input_format_type::counts:
      return load_counts(infile);
    case input_format_type::bam:
      return load_counts_BAM_se(n_threads, infile);
    default:  // case input_format_type::vals:
      return load_counts_bed_se(infile);
    }
  }();

  const auto max_observed_count = std::size(counts_hist) - 1;
  const auto n_distinct =
    std::reduce(std::cbegin(counts_hist), std::cend(counts_hist));
  const auto n_total = get_counts_from_hist(counts_hist);
  const auto gt0 = [](const auto c) { return c > 0.0; };
  const auto distinct_counts = std::ranges::count_if(counts_hist, gt0);
  if (step_size > 0)
    n_points =
      static_cast<std::uint32_t>(static_cast<std::uint64_t>(n_total) /
                                 static_cast<std::uint64_t>(step_size));
  if (verbose)
    std::print("TOTAL READS     = {}\n"
               "COUNTS_SUM      = {}\n"
               "DISTINCT READS  = {}\n"
               "DISTINCT COUNTS = {}\n"
               "MAX COUNT       = {}\n"
               "COUNTS OF 1     = {}\n"
               "N POINTS        = {}\n",
               n_reads,             //
               n_total,             //
               n_distinct,          //
               distinct_counts,     //
               max_observed_count,  //
               counts_hist[1],      //
               n_points             //
    );

  if (n_points < 2) {
    std::println("too few interpolation points; change step size");
    return EXIT_FAILURE;
  }

  const auto interps = preseq::interpolate_distinct(counts_hist, n_points);
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);
  std::println(out, "{}\t{}", "n_total", "n_distinct");
  std::println(out, "{}\t{}", 0, 0);
  std::ranges::for_each(interps, [&](const auto i) {
    std::println(out, "{}\t{}", std::round(i.first), std::round(i.second));
  });

  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
