// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#ifndef SRC_BOUND_POP_HPP_
#define SRC_BOUND_POP_HPP_

#include <span>

const auto bound_pop_about_msg =
  R"(Estimate population size (alternate method))";

auto
bound_pop_main(const std::span<char *> args) -> int;

#endif  // SRC_BOUND_POP_HPP_
