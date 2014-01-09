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

#ifndef FILE_ITERATOR_HPP
#define FILE_ITERATOR_HPP

#include <string>
#include <vector>
#include <fstream>

#include "smithlab_utils.hpp"

template <class T>
class FileIterator {
public:
  FileIterator(const std::string f, const size_t bs);
  void increment_first() {
    if (++first == buffer.end()) {
      assert(first <= last);
      refill_buffer();
    }
    assert(first <= buffer.end() && last <= buffer.end());
  }
  void increment_last() {
    assert(last < buffer.end());
    if (++last == buffer.end()) {
      assert(first <= last);
      refill_buffer();
    }
    assert(first <= buffer.end());
    assert(last <= buffer.end());
  }
  void increment() {
    increment_last();
    increment_first();
  }
  typename std::vector<T>::const_iterator get_first() const {return first;}
  typename std::vector<T>::const_iterator get_last() const {return last;}
  bool first_is_good() const {return (!in.eof() || first < buffer.end());}
  bool last_is_good() const {return (!in.eof() || last < buffer.end());}
  bool is_good() const {return first_is_good() && last_is_good();}
  
private:
  std::ifstream in;
  std::vector<T> buffer;
  typename std::vector<T>::iterator first;
  typename std::vector<T>::iterator last;
  void refill_buffer();
};

void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    std::vector<std::string> &buffer);

class GenomicRegion;
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    std::vector<GenomicRegion> &buffer);

typedef std::pair<std::string, std::string> FASTQRecord;
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    std::vector<FASTQRecord> &buffer);

class MappedRead;
void 
fill_buffer(std::ifstream &in, const size_t buffer_start, 
	    std::vector<MappedRead> &buffer);

/* THIS REFILL BUFFER IS USED WHEN INCREMENTS TO EITHER THE FIRST OR
   THE LAST CURRENTLY USED ELEMENTS IN THE BUFFER HIT THE END OF THE
   BUFFER. HOPEFULLY THE FIRST ONE WILL NOT HIT THE END BEFORE THE
   LAST: THE FIRST IS ALWAYS SUPPOSED TO BE LESS THAN OR EQUAL TO THE
   LAST.
 */
template <class T> void
FileIterator<T>::refill_buffer() {
  assert(first <= last);
  const size_t diff = last - first;
  copy(first, last, buffer.begin());
  // Not sure if the code below actualy works or is the best way to
  // grow the buffer
  // assert(diff < buffer.size());
  if (diff == buffer.size()) {
    std::vector<T> newbuff(2*buffer.size());
    copy(buffer.begin(), buffer.end(), newbuff.begin());
    buffer.swap(newbuff);
  }
  first = buffer.begin();
  last = first + diff;
  fill_buffer(in, diff, buffer);
}

template <class T>
FileIterator<T>::FileIterator(const std::string f, const size_t bs) :
  buffer(std::vector<T>(bs)) {
  in.open(f.c_str());
  if (!in) throw SMITHLABException("cannot open input file " + f);
  fill_buffer(in, 0, buffer);
  first = buffer.begin();
  last = buffer.begin();
}

#endif
