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

#include "MappedRead.hpp"
#include "smithlab_utils.hpp"

#include <cassert>
#include <fstream>
#include <tr1/unordered_map>
#include <algorithm>

using std::string;
using std::vector;

static void
find_sixth_and_seventh_whitespace(const char *buffer, 
				  size_t &sixth_ws, size_t &seventh_ws) {
  size_t ws_count = 0;
  const char *position = buffer;
  while (ws_count < 6 && *position != '\0') {
    if (*position == '\t' || *position == ' ') ++ws_count;
    if (ws_count < 6) 
      ++position;
  }
  if (ws_count != 6)
    throw SMITHLABException("malformed line (a):\n" + string(buffer));
  sixth_ws = std::distance(buffer, position);
  ++position;
  while (*position != '\0' && *position != '\t' && *position != ' ')
    ++position;
  if (*position != '\t' && *position != ' ')
    throw SMITHLABException("malformed line (b):\n" + string(buffer) + "\t\"" + *position + "\"");
  seventh_ws = std::distance(buffer, position);
}

MappedRead::MappedRead(const char *line) : r(line) {
  size_t sixth_ws = 0, seventh_ws = 0;
  find_sixth_and_seventh_whitespace(line, sixth_ws, seventh_ws);
  seq = string(line + sixth_ws + 1, seventh_ws - sixth_ws - 1);
  scr = string(line + seventh_ws + 1);
}

std::istream& 
operator>>(std::istream& the_stream, MappedRead &mr) {
  string chr, name;
  size_t start = 0ul, end = 0ul;
  char strand = '\0';
  double score;
  if (!(the_stream >> chr >> start >> end >> name >> 
	score >> strand >> mr.seq >> mr.scr))
    the_stream.setstate(std::ios::badbit);
  
  mr.r = GenomicRegion(chr, start, end, name, score, strand);
  
  char c;
  while ((c = the_stream.get()) != '\n' && the_stream);
  
  if (c != '\n')
    the_stream.setstate(std::ios::badbit);
  
  // the_stream.peek();
  if (the_stream.eof())
    the_stream.setstate(std::ios::badbit);
  
  return the_stream;
}

void
LoadMappedReadsFile(string filename, 
		    vector<MappedRead> &the_mapped_reads) {
  static const size_t buffer_size = 10000; // Magic
  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw SMITHLABException("cannot open input file " + filename);
  
  char buffer[buffer_size];
  while (!in.eof()) {
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    the_mapped_reads.push_back(MappedRead(buffer));
    in.peek();
  }
  in.close();
}

std::ostream& 
operator<<(std::ostream& the_stream, const MappedRead &mr) {
  return the_stream << mr.r << '\t' << mr.seq << '\t' << mr.scr;
}
