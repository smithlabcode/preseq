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

#include "RNG.hpp"

#include <unistd.h>
#include <cstdlib>
#include <ctime>

size_t Runif::state = 0;
bool Runif::seed_set = false;

static const size_t 
DUMMY_SEED = std::numeric_limits<size_t>::max();
static const size_t
MODULUS_MASK = static_cast<size_t>(-1);
static const double
DOUBLE_DENOMINATOR = static_cast<double>(std::numeric_limits<int>::max());

#include <iostream>

static inline size_t
rng_integer(size_t &state) {
  // TOO MUCH MAGIC:
  // state = (1103515245*state + 12345) & MODULUS_MASK;
  state = rand();
  return state;
}

static double
rng_double(size_t &state) {
  return rng_integer(state)/DOUBLE_DENOMINATOR;
}

Runif::Runif(size_t seed) {
  if (seed != DUMMY_SEED)
    state = seed;
  else if (!seed_set)
    state = time(0) + getpid();
  srand(state);
}

int
Runif::runif(int min_val, int max_val) const {
  return min_val + (rng_integer(state) % (max_val - min_val));
}

size_t
Runif::runif(size_t min_val, size_t max_val) const {
  return min_val + (rng_integer(state) % (max_val - min_val));
}

double 
Runif::runif(double min_val, double max_val) const {
  return min_val + rng_double(state)*(max_val - min_val);
}
