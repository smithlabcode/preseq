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

#include "bound_pop.hpp"

#include "common.hpp"
#include "lnfact.hpp"
#include "load_data_for_complexity.hpp"
#include "moment_sequence.hpp"

#include "CLI11/CLI11.hpp"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>  // IWYU pragma: keep
#include <numeric>
#include <print>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

// bounding n_0
auto
bound_pop_main(int argc, char *argv[]) -> int {  // NOLINT (*-avoid-c-arrays)
  try {
    const auto normalize =
      [](auto &x) {  // cppcheck-suppress[constParameterReference]
        const auto d = std::reduce(std::cbegin(x), std::cend(x));
        std::ranges::transform(x, std::begin(x),
                               [&](const auto y) { return y / d; });
      };

    bool verbose{false};
    bool paired_end{false};
    bool quick_mode{false};

    std::string infile;
    std::string outfile;
    std::string histogram_outfile;

    // NOLINTBEGIN(*-avoid-magic-numbers)
    std::size_t max_num_points = 10;
    double tolerance = 1e-20;
    std::size_t n_bootstraps = 500;
    double c_level = 0.95;
    std::size_t max_iter = 100;
    std::uint32_t seed = 408;
    // NOLINTEND(*-avoid-magic-numbers)

#ifdef HAVE_HTSLIB
    std::uint32_t n_threads{1};
#endif

    CLI::App app{rlstrip(bound_pop_about_msg)};
    argv = app.ensure_utf8(argv);
    app.formatter(std::make_shared<preseq_formatter>());
    app.usage("\nUsage: preseq bound_pop [OPTIONS]");
    // if (argc >= 3)
    //   app.footer(rlstrip(description));

    // clang-format off
    app.add_option("INPUT", infile, "input file name")
      ->required()
      ->option_text(" ")
      ->check(CLI::ExistingFile)
      // ->check(CLI::ReadPermission)
      ;
    app.add_option("-o,--output", outfile, "output file");
    app.add_option("-m,--max-points", max_num_points,
                   "maximum number of points in quadrature estimates");
    app.add_option("-t,--tolerance", tolerance, "numerical tolerance");
    app.add_option("-n,--bootstraps", n_bootstraps, "number of bootstraps");
    app.add_option("-c,--ci-level", c_level, "level for confidence intervals");
    app.add_option("-r,--seed", seed, "seed for random number generator");
    app.add_flag("-p,--paired-end", paired_end, "input is paired end read file");
    app.add_flag("-q,--quick", quick_mode, "no bootstraps when making estimates");
    app.add_flag("-v,--verbose", verbose, "print moments and boostraps with output");
    // clang-format on

    if (argc < 3) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const auto input_format = get_input_format_type(infile);
    if (is_unknown(input_format)) {
      std::println("unknown input format");
      return EXIT_FAILURE;
    }

    const auto [n_obs, counts_hist] = [&] {
      if (is_hist(input_format))
        return load_histogram(infile);
      if (is_counts(input_format))
        return load_counts(infile);
#ifdef HAVE_HTSLIB
      if (is_bam(input_format))
        return paired_end ? load_counts_BAM_pe(n_threads, infile)
                          : load_counts_BAM_se(n_threads, infile);
#endif
      //  if (is_bed(input_format))
      return paired_end ? load_counts_bed_pe(infile)
                        : load_counts_bed_se(infile);
    }();

    const double distinct_obs =
      std::accumulate(std::cbegin(counts_hist), std::cend(counts_hist), 0.0);

    std::vector<double> measure_moments;
    // mu_r = (r + 1)! n_{r+1} / n_1
    for (auto i = 1u; i < std::size(counts_hist) && counts_hist[i]; ++i) {
      const auto mm = std::exp(lnfact(i + 1) + std::log(counts_hist[i]) -
                               std::log(counts_hist[1]));
      if (!std::isfinite(mm))
        break;
      measure_moments.push_back(mm);
    }

    if (!histogram_outfile.empty())
      report_histogram(histogram_outfile, counts_hist);

    std::vector<nlohmann::json> bootstraps;
    nlohmann::json output;

    if (quick_mode) {
      if (std::size(measure_moments) > 2 * max_num_points)
        measure_moments.resize(2 * max_num_points);

      auto n_points = ensure_pos_def_mom_seq(measure_moments, tolerance);
      MomentSequence obs_mom_seq(measure_moments);

      std::vector<double> points, weights;
      obs_mom_seq.lower_quadrature_rules(n_points, tolerance, max_iter, points,
                                         weights);
      normalize(weights);

      const auto n_1 = counts_hist[1];
      const auto term = [n_1](const auto w, const auto p) {
        return n_1 * w / p;
      };
      auto estimated_unobs =
        std::inner_product(std::cbegin(weights), std::cend(weights),
                           std::cbegin(points), 0.0, std::plus<>(), term);
      estimated_unobs = std::max(estimated_unobs, 0.0) + distinct_obs;
      if (estimated_unobs == distinct_obs)
        n_points = 0;

      output = {
        {"quadrature_estimated_unobs", estimated_unobs},
        {"n_points", n_points},
      };
    }
    else {
      // do bootstraps
      std::vector<double> quad_estimates;

      std::mt19937 rng(seed);  // setup rng

      // hist may be sparse, to speed up bootstrapping
      // sample only from positive entries
      std::vector<std::size_t> counts_hist_distinct_counts;
      std::vector<double> distinct_counts_hist;
      for (std::size_t i = 0; i < std::size(counts_hist); ++i)
        if (counts_hist[i] > 0) {
          counts_hist_distinct_counts.push_back(i);
          distinct_counts_hist.push_back(counts_hist[i]);
        }

      for (auto i = 0u;
           i < max_iter && std::size(quad_estimates) < n_bootstraps; ++i) {
        std::vector<double> sample_hist;
        resample_hist(rng, counts_hist_distinct_counts, distinct_counts_hist,
                      sample_hist);
        const double sampled_distinct = std::accumulate(
          std::cbegin(sample_hist), std::cend(sample_hist), 0.0);

        // initialize moments, 0-th moment is 1
        std::vector<double> bootstrap_moments(1, 1.0);
        // moments[r] = (r + 1)! n_{r+1} / n_1
        for (std::size_t j = 0; j < 2 * max_num_points; ++j)
          bootstrap_moments.push_back(std::exp(lnfact(j + 3) +
                                               std::log(sample_hist[j + 2]) -
                                               std::log(sample_hist[1])));
        const auto n_points = std::min(
          ensure_pos_def_mom_seq(bootstrap_moments, tolerance), max_num_points);

        MomentSequence bootstrap_mom_seq(bootstrap_moments);

        std::vector<double> points;
        std::vector<double> weights;
        bootstrap_mom_seq.lower_quadrature_rules(n_points, tolerance, max_iter,
                                                 points, weights);
        normalize(weights);

        const auto n_1 = counts_hist[1];
        const auto term = [n_1](const auto w, const auto p) {
          return n_1 * w / p;
        };
        auto estimated_unobs =
          std::inner_product(std::cbegin(weights), std::cend(weights),
                             std::cbegin(points), 0.0, std::plus{}, term);
        estimated_unobs = std::max(estimated_unobs, 0.0) + sampled_distinct;

        if (verbose)
          bootstraps.push_back(nlohmann::json({
            {"bootstrapped_moments", bootstrap_moments},
            {"alpha", bootstrap_mom_seq.alpha},
            {"beta", bootstrap_mom_seq.alpha},
            {"points", points},
            {"weights", weights},
            {"estimated_unobs", estimated_unobs},
          }));

        quad_estimates.push_back(estimated_unobs);
      }

      double median_estimate{};
      double lower_ci{};
      double upper_ci{};
      median_and_ci(quad_estimates, c_level, median_estimate, lower_ci,
                    upper_ci);
      output = {
        {"median_estimated_unobs", median_estimate},
        {"lower_ci", lower_ci},
        {"upper_ci", upper_ci},
      };
    }

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    if (!outfile.empty() && !out)
      throw std::runtime_error("failed to open output file: " + outfile);

    output["total_observations"] = n_obs;
    output["distinct_observations"] = distinct_obs;
    output["max_count"] = std::size(counts_hist) - 1;
    output["observed_moments"] = measure_moments;
    if (verbose && !quick_mode)
      output["bootstraps"] = bootstraps;
    std::println(out, "{}", output.dump(4));
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
