/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
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

#ifndef SIM_UTILS_HPP
#define SIM_UTILS_HPP

#include <string>
#include <vector>
#include "RNG.hpp"

void
sequence_to_consensus_matrix(const std::string &sequence, 
			     std::vector<std::vector<double> > &matrix);

void
add_sequencing_errors(const Runif &rng, const double n_errors,
		      std::vector<std::vector<double> > &matrix);

void
add_sequencing_errors(const double n_errors,
		      std::vector<std::vector<double> > &matrix);

void
call_bases_solexa(const std::vector<std::vector<double> > &matrix,
		  std::string &sequence);

#endif
