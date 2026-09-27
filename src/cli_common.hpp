// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_CLI_COMMON_HPP_
#define SRC_CLI_COMMON_HPP_

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
public:
  auto
  make_option_desc(const CLI::Option *opt) const -> std::string override {
    static constexpr auto max_descr_width = 50;
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
};

#endif  // SRC_CLI_COMMON_HPP_
