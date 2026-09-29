// SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

#ifndef SRC_GC_EXTRAP_HPP_
#define SRC_GC_EXTRAP_HPP_

#include <span>

static constexpr auto gc_extrap_about_msg =
  R"(Estimate the fraction of the genome covered by reads)";

static constexpr auto gc_extrap_footer_msg =
  R"(This approach is described in Daley & Smith (2014). The method is the same as
for lc_extrap: using rational function approximation to a power-series
expansion for the number of "unobserved" bases in the initial sample. The
gc_extrap method is adapted to deal with individual nucleotides rather than
distinct reads.)";

auto
gc_extrap_main(const std::span<char *> args) -> int;

#endif  // SRC_GC_EXTRAP_HPP_
