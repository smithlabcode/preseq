// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_CLI_COMMON_HPP_
#define SRC_CLI_COMMON_HPP_

#define CLI11_ENABLE_EXTRA_VALIDATORS 1

#include <CLI11/CLI11.hpp>

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

class preseq_formatter : public CLI::Formatter {
  static constexpr auto total_output_width_ = 80;

public:
  preseq_formatter() : Formatter() {
    CLI::FormatterBase::enable_default_flag_values_ = false;
  }
  auto
  make_option_desc(const CLI::Option *opt) const -> std::string override {
    const auto max_descr_width = right_column_width_;
    std::istringstream iss{opt->get_description()};
    const std::vector<std::string> words{
      std::istream_iterator<std::string>{iss}, {}};
    std::string r{words[0]};
    std::uint32_t width = std::size(words[0]);
    for (auto i = 1u; i < std::size(words); ++i) {
      if (width == 0 || width + std::size(words[i]) < max_descr_width) {
        r += ' ';
        ++width;
      }
      else {
        r += '\n';
        width = 0;
      }
      r += words[i];
      width += std::size(words[i]);
    }
    return r;
  }
  auto
  column_width(const std::size_t val) -> void {
    column_width_ = val;
    right_column_width_ = total_output_width_ - column_width_;
  }
};

#endif  // SRC_CLI_COMMON_HPP_
