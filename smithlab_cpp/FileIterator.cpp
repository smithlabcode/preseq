/*
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

#include "FileIterator.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include <iostream>

using std::istream;
using std::vector;
using std::string;

/*  THIS FUNCTION FILLS A BUFFER FOR GenomicRegion OBJECTS
 */
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    vector<GenomicRegion> &buffer) {
  GenomicRegion tmp;
  size_t i = buffer_start;
  assert(buffer_start <= buffer.size());
  for (; i != buffer.size() && !in.eof(); ++i) {
    in >> tmp;
    buffer[i].swap(tmp);
    in.peek();
  }
  if (i < buffer.size())
    buffer.erase(buffer.begin() + i, buffer.end());
}

/*  THIS FUNCTION FILLS A BUFFER FOR THE ACTUAL READS, REPRESENTED AS
 *  STRINGS, AND MUST BE IN A FASTA FORMAT FILE
 */
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    vector<string> &buffer) {
  string tmp;
  size_t i = buffer_start;
  for (; i != buffer.size() && !in.eof(); ++i) {
    // the read name...
    in >> tmp; // DANGER: assumes that the name of the read has no
	       // spaces in it!!
    // the read itself:
    in >> buffer[i];
    in.peek();
  }
  if (i < buffer.size())
    buffer.erase(buffer.begin() + i, buffer.end());
}


/*  THIS FUNCTION FILLS A BUFFER FOR THE ACTUAL READS, REPRESENTED AS
 *  RECORDS IN A FASTQ FILE, INCLUDING THE QUALITY SCORES
 */
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    vector<FASTQRecord> &buffer) {
  string tmp, read_seq, scores_seq;
  size_t i = buffer_start;
  for (; i != buffer.size() && !in.eof(); ++i) {
    // the read name...
    in >> tmp; // DANGER: assumes that the name of the read has no
	       // spaces in it!!
    // the read itself:
    in >> read_seq;
    in >> tmp;
    in >> scores_seq;
    buffer[i] = std::make_pair(read_seq, scores_seq);
    in.peek();
  }
  if (i < buffer.size())
    buffer.erase(buffer.begin() + i, buffer.end());
}


/*  THIS FUNCTION FILLS A BUFFER FOR THE ACTUAL READS, REPRESENTED AS
 *  RECORDS IN A FASTQ FILE, INCLUDING THE QUALITY SCORES
 */
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    vector<MappedRead> &buffer) {
  size_t i = buffer_start;
  for (; i != buffer.size() && !in.eof(); ++i) {
    in >> buffer[i];
    in.peek();
  }
  if (i < buffer.size())
    buffer.erase(buffer.begin() + i, buffer.end());
}
