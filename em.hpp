#ifndef EM_HPP
#define EM_HPP

#include "rmap_utils.hpp"

void
EM(const std::vector<double> &values,
   const size_t max_iterations, const double tolerance,
   bool VERBOSE, double &fg_distro, double &mid_distro,
   double &bg_distro, std::vector<double> &mixing);

#endif
