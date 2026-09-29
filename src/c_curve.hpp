// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#ifndef SRC_C_CURVE_HPP_
#define SRC_C_CURVE_HPP_

#include <span>

static constexpr auto c_curve_about_msg =
  R"(Generate the interpolated complexity curve (no estimation))";

auto
c_curve_main(const std::span<char *> args) -> int;

#endif  // SRC_C_CURVE_HPP_
