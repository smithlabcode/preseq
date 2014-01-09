/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RNG_HPP
#define RNG_HPP

#include <limits>
#include <cstdlib>

class Runif {
public:
  Runif(size_t seed = std::numeric_limits<size_t>::max());
  ~Runif() {}
  int runif(int min_val, int max_val) const;
  size_t runif(size_t min_val, size_t max_val) const;
  double runif(double min_val, double max_val) const;  
private:
  static size_t state;
  static bool seed_set;
};

#endif
