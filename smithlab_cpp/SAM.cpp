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

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "SAM.hpp"
#include "smithlab_utils.hpp"
#include "MappedRead.hpp"
#include "sam.h"

using std::string;
using std::vector;
using std::cerr;
using std::endl;

SAMReader::SAMReader(const string fn, const string mapper_used) :
  filename(fn), mapper(mapper_used), file_handler(NULL),
  algn_p(NULL), GOOD(false) 
{
  if (mapper != "bsmap" && mapper != "bismark" 
      && mapper != "bs_seeker" && mapper != "general")
    throw SMITHLABException("Mapper unsupported:" + mapper);

  const string ext_name = filename.substr(filename.find_last_of('.'));
  mode = ext_name == ".bam" ? "rb" : "r";
  
  if ((file_handler = samopen(filename.c_str(), mode.c_str(), NULL))
      == NULL)
  {
    cerr << "Fail to open SAM/BAM file " << filename << endl;
    exit(-1);
  }
  
  algn_p = bam_init1();
  GOOD = true;
}

SAMReader::~SAMReader()
{
  close();
}

void
SAMReader::close()
{
  if (algn_p)
  {
    bam_destroy1(algn_p);
    algn_p = NULL;
  }
  if (file_handler)
  {
    samclose(file_handler);
    file_handler = NULL;
    filename = "";
    mode = "";
    GOOD = false;
  }
}

bool
SAMReader::get_SAMRecord(const string &str, SAMRecord &samr)
{
  if (mapper == "bsmap")
    return get_SAMRecord_bsmap(str, samr);
  else if (mapper == "bismark")
    return get_SAMRecord_bismark(str, samr);
  else if (mapper == "bs_seeker")
    return get_SAMRecord_bsseeker(str, samr);
  else if (mapper == "general")
    return get_SAMRecord_general(str, samr);
  else
    GOOD = false;
  return false;
}

SAMReader&
operator>>(SAMReader &sam_stream, SAMRecord& samr)
{
  if (samread(sam_stream.file_handler, sam_stream.algn_p) >= 0)
  {
    char *str = bam_format1_core(sam_stream.file_handler->header,
                                       sam_stream.algn_p,
                                       sam_stream.file_handler->type>>2&3);
    sam_stream.GOOD = sam_stream.get_SAMRecord(str, samr);
    free(str);
  }
  else
    sam_stream.GOOD = false;

  return sam_stream;
}

/////////////////////////////////////////////
//// general facility for SAM format
/////////////////////////////////////////////

void static
apply_CIGAR(const string &seq, const string &qual,
            const string &CIGAR, string &new_seq, string &new_qual)
{
    assert(seq.size() == qual.size());
    assert(new_seq.size() == 0 && new_qual.size() == 0);
    size_t n;
    char op;
    size_t i = 0;

    std::istringstream iss(CIGAR);
    while (iss >> n >> op)
    {
        switch (op)
        {
        case 'M':
            new_seq += seq.substr(i, n);
            new_qual += qual.substr(i, n);
            i += n;
            break;
        case 'I':
            i += n;
            break;
        case 'D':
            new_seq += string(n, 'N');
            new_qual += string(n, 'B');
            break;
        case 'S':
            i += n;
            break;
        case 'H':
            ;
            break;
        case 'P':
            ;
            break;
        case '=':
	    new_seq += seq.substr(i, n);
            new_qual += qual.substr(i, n);
            i += n;
            break;
        case 'X':
            new_seq += seq.substr(i, n);
            new_qual += qual.substr(i, n);
            i += n;
            break;
        }
    }
    // Sum of lengths of the M/I/S/=/X operations 
    // shall equal the length of seq.

    assert(i == seq.length());
    assert(new_seq.size() == new_qual.size());
}

class FLAG 
{
public:  
  FLAG(const size_t f) : flag(f) {}
  bool is_pairend() const {return flag & 0x1;}
  bool is_singlend() const {return !(is_pairend());}
  bool is_mapping_paired() const {return flag & 0x2;}
  bool is_unmapped() const {return flag & 0x4;}
  bool is_mapped() const {return !(is_unmapped());}
  bool is_revcomp() const {return flag & 0x10;}
  bool is_Trich() const {return flag & 0x40;}
  bool is_Arich() const {return flag & 0x80;}
  bool is_secondary() const {return flag & 0x100;}
  bool is_primary() const {return !(is_secondary());}

private:
  size_t flag;
};

////////////////////////////////////////
// BSMAP
////////////////////////////////////////

class BSMAPFLAG : public FLAG
{
public:  
  BSMAPFLAG(const size_t f) : FLAG(f) {}
};

inline static void
bsmap_get_strand(const string &strand_str, string &strand, string &bs_forward)
{
    strand = strand_str.substr(5, 1);
    bs_forward = strand_str.substr(6, 1);
    if (bs_forward == "-") strand = strand == "+" ? "-" : "+";
}

bool
SAMReader::get_SAMRecord_bsmap(const string &str, SAMRecord &samr)
{
/////
  cerr << "WARNING: "<< "[BSMAP Converter] test version: may contain bugs" << endl;
/////

  string name, chrom, CIGAR, mate_name, seq, qual, strand_str, mismatch_str;
  size_t flag, start, mapq_score, mate_start;
  int seg_len;

  std::istringstream iss(str);
  if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual
        >> mismatch_str >> strand_str))
  {
    GOOD = false;
    throw SMITHLABException("malformed line in bsmap SAM format:\n" + str);
  }
  
  BSMAPFLAG Flag(flag);

  samr.mr.r.set_chrom(chrom);
  samr.mr.r.set_start(start - 1);
  samr.mr.r.set_name(name);
  samr.mr.r.set_score(atoi(mismatch_str.substr(5).c_str()));
  
  string strand, bs_forward;
  bsmap_get_strand(strand_str, strand, bs_forward);
  samr.mr.r.set_strand(strand[0]);

  string new_seq, new_qual;
  apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);

  samr.mr.r.set_end(samr.mr.r.get_start() + new_seq.size()); 
  samr.mr.seq = new_seq;
  samr.mr.scr = new_qual;

  samr.is_Trich = Flag.is_Trich();
  samr.is_mapping_paired = Flag.is_mapping_paired();
  return GOOD;
}

////////////////////////////////////////
// Bismark
////////////////////////////////////////

class BISMARKFLAG : public FLAG
{
public:  
  BISMARKFLAG(const size_t f) : FLAG(f) {}
  bool is_Trich() const {return is_pairend() ? FLAG::is_Trich() : true;}
};

static size_t
get_mismatch_bismark(const string &edit_distance_str,
                     const string &meth_call_str)
{
  /*
  the result of this function might not be accurate, because if a sequencing
  error occurs on a cytosine, then it probably will be reported as a convertion
  */
  size_t edit_distance;
  edit_distance = atoi(edit_distance_str.substr(5).c_str());

  int convert_count = 0;
  const char *temp = meth_call_str.substr(5).c_str();
  while(*temp != '\0') {
    if (*temp == 'x' || *temp == 'h' || *temp == 'z')
      ++convert_count;
    ++temp;
  }

  return edit_distance - convert_count;
}

bool
SAMReader::get_SAMRecord_bismark(const string &str, SAMRecord &samr)
{
  string name, chrom, CIGAR, mate_name, seq, qual, strand_str,
    edit_distance_str, mismatch_str, meth_call_str,
    read_conv_str, genome_conv_str;
  size_t flag, start, mapq_score, mate_start;
  int seg_len;
  
  std::istringstream iss(str);
  if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual
        >> edit_distance_str >> mismatch_str >> meth_call_str
        >> read_conv_str >> genome_conv_str))
  {
    GOOD = false;
    throw SMITHLABException("malformed line in bismark SAM format:\n" + str);
  }

  BISMARKFLAG Flag(flag);

  samr.mr.r.set_chrom(chrom);
  samr.mr.r.set_start(start - 1);
  samr.mr.r.set_name(name);
  samr.mr.r.set_score(get_mismatch_bismark(edit_distance_str, meth_call_str));
  samr.mr.r.set_strand(Flag.is_revcomp() ? '-' : '+');
  
  string new_seq, new_qual;
  apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);

  if (Flag.is_revcomp())
  {
    revcomp_inplace(new_seq);
    std::reverse(new_qual.begin(), new_qual.end());
  }

  samr.mr.r.set_end(samr.mr.r.get_start() + new_seq.size()); 
  samr.mr.seq = new_seq;
  samr.mr.scr = new_qual;

  samr.is_Trich = Flag.is_Trich();
  samr.is_mapping_paired = Flag.is_mapping_paired();

// /////
//   cerr << "check 1: "<< (samr.is_Trich ? "T-rich" : "A-rich") << "\t"
//        << samr.is_mapping_paired << endl
//        << samr.mr << endl;
// /////
  return GOOD;
}

////////////////////////////////////////
// BS Seeker
////////////////////////////////////////

class BSSEEKERFLAG : public FLAG
{
public:  
  BSSEEKERFLAG(const size_t f) : FLAG(f) {}
  // pair-end:
  //  if T-rich mate is +, then both mates are +;
  //  if T-rich mate is -, then both mates are -;
  // single-end:
  //  0 for +; 16 for -.
  bool is_Trich() const {
    return FLAG::is_pairend() ? FLAG::is_Trich() : true;
  }
  bool is_Arich() const {
    return FLAG::is_pairend() && FLAG::is_Arich();
  }
  bool is_revcomp() const {
    if (FLAG::is_pairend())
      return FLAG::is_revcomp() ? is_Trich() : is_Arich();
    else
      return FLAG::is_revcomp();
  }

};

bool
SAMReader::get_SAMRecord_bsseeker(const string &str, SAMRecord &samr)
{
  
  string name, chrom, CIGAR, mate_name, seq, qual, orientation_str,
      conversion_str, mismatch_str, mismatch_type_str, seq_genome_str;
  size_t flag, start, mapq_score, mate_start;
  int seg_len;
  
  std::istringstream iss(str);
  if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual
        >> orientation_str >> conversion_str >> mismatch_str
        >> mismatch_type_str >> seq_genome_str))
  {
    GOOD = false;
    throw SMITHLABException("malformed line in bs_seeker SAM format:\n" + str);
  }

  // bs_seeker also doesn't keep sequencing quality information?
  qual = string(seq.size(), 'h');

  BSSEEKERFLAG Flag(flag);
  
  samr.mr.r.set_chrom(chrom);
  samr.mr.r.set_start(start - 1);
  samr.mr.r.set_name(name);
  samr.mr.r.set_score(atoi(mismatch_str.substr(5).c_str())); 
  samr.mr.r.set_strand(Flag.is_revcomp() ? '-' : '+'); 

  string new_seq, new_qual;
  apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);
  
  if (Flag.is_revcomp())
  {
    revcomp_inplace(new_seq);
    std::reverse(new_qual.begin(), new_qual.end());
  }
  
  samr.mr.r.set_end(samr.mr.r.get_start() + new_seq.size()); 
  samr.mr.seq = new_seq;
  samr.mr.scr = new_qual;

  samr.is_Trich = Flag.is_Trich();
  samr.is_mapping_paired = Flag.is_mapping_paired();
  
  return GOOD;
}

////////////////////////////////////////
// general: for non-BSseq mappers
////////////////////////////////////////

class GENERALFLAG : public FLAG
{
public:  
  GENERALFLAG(const size_t f) : FLAG(f) {}
  bool is_Trich() const {return is_pairend() ? FLAG::is_Trich() : true;}
};

bool
SAMReader::get_SAMRecord_general(const string &str, SAMRecord &samr)
{
  string name, chrom, CIGAR, mate_name, seq, qual;
  size_t flag, start, mapq_score, mate_start;
  int seg_len;
  
  std::istringstream iss(str);
  if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual))
  {
    GOOD = false;
    throw SMITHLABException("malformed line in SAM format:\n" + str);
  }

 
  GENERALFLAG Flag(flag);
  samr.is_primary = Flag.is_primary();
  samr.is_mapped = Flag.is_mapped();

  samr.seg_len = seg_len;
  samr.mr.r.set_name(name);

  if(samr.is_primary && samr.is_mapped){
    samr.mr.r.set_chrom(chrom);
    samr.mr.r.set_start(start - 1);
    samr.mr.r.set_score(0);
    samr.mr.r.set_strand(Flag.is_revcomp() ? '-' : '+');
  
    string new_seq, new_qual;
    apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);

    
    if (Flag.is_revcomp())
    {
      revcomp_inplace(new_seq);
      std::reverse(new_qual.begin(), new_qual.end());
    }
    

    samr.mr.r.set_end(samr.mr.r.get_start() + new_seq.size());
     
    samr.mr.seq = new_seq;
    samr.mr.scr = new_qual;

  }
  samr.is_Trich = Flag.is_Trich();
  samr.is_mapping_paired = Flag.is_mapping_paired();

  return GOOD;
}
