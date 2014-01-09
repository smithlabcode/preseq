/*
 *    Part of SMITHLAB_CPP software
 *
 *    Copyright (C) 2010 University of Southern California and
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

#ifndef MAPPED_READ_HPP
#define MAPPED_READ_HPP

#include "GenomicRegion.hpp"

struct MappedRead {
  MappedRead() {}
  MappedRead(const char *line);
  GenomicRegion r;
  std::string seq;
  std::string scr;
};

void
LoadMappedReadsFile(std::string filename, 
		    std::vector<MappedRead> &the_mapped_reads);

std::istream& 
operator>>(std::istream& the_stream, MappedRead &mr);
std::ostream& 
operator<<(std::ostream& the_stream, const MappedRead &mr);

#endif
