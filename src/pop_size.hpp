// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#ifndef SRC_POP_SIZE_HPP_
#define SRC_POP_SIZE_HPP_

#include <span>

static constexpr auto pop_size_about_msg =
  R"(Estimate population size (preferred method))";

static constexpr auto pop_size_footer_msg =
  R"(Estimate the total population size using the approach described in Daley & Smith
(2013), extrapolating to very long range. Default parameters assume that the
initial sample represents at least 1e-9 of the population, which is sufficient
for every example application we have seen.)";

auto
pop_size_main(const std::span<char *> args) -> int;

#endif  // SRC_POP_SIZE_HPP_
