/*
 *    Part of SMITHLAB_CPP software
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Meng Zhou, Qiang Song, Andrew Smith
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

#ifndef SAM_HPP
#define SAM_HPP

#include "smithlab_utils.hpp"
#include "MappedRead.hpp"
#include "sam.h"

#include <string>
#include <vector>
#include <fstream>

struct SAMRecord {
  MappedRead mr;
  bool is_Trich;
  bool is_mapping_paired;
  bool is_primary;
  bool is_mapped;
  int seg_len;
};

class SAMReader
{
public:  
  SAMReader(const std::string filename, const std::string mapper);
  ~SAMReader();
  
  friend SAMReader& 
  operator>>(SAMReader& sam_stream, SAMRecord &samr);
  
  const bool is_good() const {return GOOD;}

private:  
  void
  close();

  // internal methods
  bool
  get_SAMRecord(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_bsmap(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_bismark(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_bsseeker(const std::string&, SAMRecord&);
  bool
  get_SAMRecord_general(const std::string&, SAMRecord&);
  
  // data
  std::string filename;
  std::string mapper;
  std::string mode;

  samfile_t*  file_handler;
  bam1_t *algn_p;
  bool GOOD;
};

#endif
