// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#ifndef SRC_LC_EXTRAP_HPP_
#define SRC_LC_EXTRAP_HPP_

#include <span>

constexpr auto lc_extrap_about_msg =
  R"(Estimate an extrapolated complexity curve)";

constexpr auto lc_extrap_footer_msg =
  R"(This is the approach described in Daley & Smith (2013). The method applies
rational function approximation via continued fractions with the original goal
of estimating the number of distinct reads that a sequencing library would yield
upon deeper sequencing.  This method has been used for many different purposes
since then.
)";

auto
lc_extrap_main(const std::span<char *> args) -> int;

#endif  // SRC_LC_EXTRAP_HPP_
