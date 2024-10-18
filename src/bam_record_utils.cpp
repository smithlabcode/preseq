/* Copyright (C) 2020-2023 Masaru Nakajima and Andrew D. Smith
 *
 * Authors: Masaru Nakajima and Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "bam_record_utils.hpp"

#include <htslib/sam.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>

#include "dnmt_error.hpp"
#include "smithlab_utils.hpp"

using std::runtime_error;
using std::string;
using std::stringstream;
using std::to_string;
using std::vector;

using bamxx::bam_header;
using bamxx::bam_rec;

/// functions in place of undefd macro
static inline bool
bam_is_rev(const bam1_t *b) {
  return (b->core.flag & BAM_FREVERSE) != 0;
}

static inline char *
bam_get_qname(const bam1_t *b) {
  return reinterpret_cast<char *>(b->data);
}

static inline uint32_t *
bam_get_cigar(const bam1_t *b) {
  return reinterpret_cast<uint32_t *>(b->data + b->core.l_qname);
}

static inline uint8_t *
bam_get_seq(const bam1_t *b) {
  // start of data +  bytes for query/read name+ bytes for cigar
  return b->data + b->core.l_qname + (b->core.n_cigar << 2);
}

static inline uint8_t *
bam_get_qual(const bam1_t *b) {
  return b->data +                     // start of data
         b->core.l_qname +             // bytes for query name
         (b->core.n_cigar << 2) +      // bytes for cigar
         ((b->core.l_qseq + 1) >> 1);  // bytes for packed query/read
}

static inline uint8_t *
bam_get_aux(const bam1_t *b) {
  return b->data + b->core.l_qname + (b->core.n_cigar << 2) +
         ((b->core.l_qseq + 1) >> 1) + b->core.l_qseq;
}

// static inline int
// bam_get_l_aux(const bam1_t *b) {
//   return b->l_data - (b->core.l_qname + (b->core.n_cigar << 2) +
//                       ((b->core.l_qseq + 1) >> 1) + b->core.l_qseq);
// }

// static inline bool
// bam_same_orientation(const bam1_t *a, const bam1_t *b) {
//   return ((a->core.flag ^ b->core.flag) & BAM_FREVERSE) != 0;
// }

static void
roundup_to_power_of_2(uint32_t &x) {
  bool k_high_bit_set = (x >> (sizeof(uint32_t) * 8 - 1)) & 1;
  if (x > 0) {
    uint8_t size = sizeof(uint32_t);
    --x;
    x |= x >> (size / 4);
    x |= x >> (size / 2);
    x |= x >> (size);
    x |= x >> (size * 2);
    x |= x >> (size * 4);
    x += !k_high_bit_set;
  }
  else
    x = 0;
}

static int
sam_realloc_bam_data(bam1_t *b, size_t desired) {
  uint32_t new_m_data = desired;
  roundup_to_power_of_2(new_m_data);
  if (new_m_data < desired) {
    errno = ENOMEM;  // (from sam.c) Not strictly true but we can't
                     // store the size
    return -1;
  }
  uint8_t *new_data = nullptr;
  if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
    new_data = static_cast<uint8_t *>(realloc(b->data, new_m_data));
  }
  else {
    new_data = static_cast<uint8_t *>(malloc(new_m_data));
    if (new_data != nullptr) {
      if (b->l_data > 0)
        std::copy_n(b->data,
                    (static_cast<uint32_t>(b->l_data) < b->m_data) ? b->l_data
                                                                   : b->m_data,
                    new_data);
      bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
    }
  }
  if (!new_data) return -1;
  b->data = new_data;
  b->m_data = new_m_data;
  return 0;
}

// static int
// sam_realloc_bam_data(bam1_t *b, size_t desired) {
//   /* returns flag: either 0 for success or -1 for error (unable to
//      allocate desired memory) */
//   uint32_t new_m_data = desired;
//   roundup_to_power_of_2(new_m_data);
//   if (new_m_data < desired)  return -1;
//   uint8_t *new_data = (uint8_t *)realloc(b->data, new_m_data);
//   if (!new_data) return -1;
//   // ADS: what would be the state of members below if -1 was returned?
//   b->data = new_data;
//   b->m_data = new_m_data;
//   return 0;
// }

// static inline void
// bam_copy_core(const bam1_t *a, bam1_t *b) {
//   // ADS: prepared for a possibly more efficient block copy to assign
//   // all variables at once, like this: b->core = a->core;
//   b->core.pos = a->core.pos;
//   b->core.tid = a->core.tid;
//   b->core.bin = a->core.bin;
//   b->core.qual = a->core.qual;
//   b->core.l_extranul = a->core.l_extranul;
//   b->core.flag = a->core.flag;
//   b->core.l_qname = a->core.l_qname;
//   b->core.n_cigar = a->core.n_cigar;
//   b->core.l_qseq = a->core.l_qseq;
//   b->core.mtid = a->core.mtid;
//   b->core.mpos = a->core.mpos;
//   b->core.isize = a->core.isize;
// }

static inline void
bam_set1_core(bam1_core_t &core, const size_t l_qname, const uint16_t flag,
              const int32_t tid, const hts_pos_t pos, const uint8_t mapq,
              const size_t n_cigar, const int32_t mtid, const hts_pos_t mpos,
              const hts_pos_t isize, const size_t l_seq,
              const size_t qname_nuls) {
  /* ADS: These are used in `hts_reg2bin` from `htslib/hts.h` and
     likely mean "region to bin" for indexing */
  /* MN: hts_reg2bin categorizes the size of the reference region.
     Here, we use the numbers used in htslib/cram/cram_samtools.h */
  static const int min_shift = 14;
  static const int n_lvls = 5;

  core.pos = pos;
  core.tid = tid;
  core.bin = hts_reg2bin(pos, pos + isize, min_shift, n_lvls);
  core.qual = mapq;
  core.l_extranul = qname_nuls - 1;
  core.flag = flag;
  core.l_qname = l_qname + qname_nuls;
  core.n_cigar = n_cigar;
  core.l_qseq = l_seq;
  core.mtid = mtid;
  core.mpos = mpos;
  core.isize = isize;
}

static inline int
bam_set1_wrapper(bam1_t *bam, const size_t l_qname, const char *qname,
                 const uint16_t flag, const int32_t tid, const hts_pos_t pos,
                 const uint8_t mapq, const size_t n_cigar,
                 const uint32_t *cigar, const int32_t mtid,
                 const hts_pos_t mpos, const hts_pos_t isize,
                 const size_t l_seq, const size_t l_aux) {
  /* This is based on how assignment is done in the `bam_set1`
     function defined in `sam.c` from htslib */

  /* This modification assigns variables of bam1_t struct but not the
   * query/read sequence.
   *
   * many checks in bam_set1 have been removed because they are
   * checked in code that calls this function, mostly because the
   * values come from valid `bam1_t` instances, so have already been
   * validated.
   *
   * Assumptions:
   * cigar: correct format and matching length
   * rlen = isize
   * qlen = l_seq
   * l_qname <= 254
   * HTS_POS_MAX - rlen > pos
   * Where HTS_POS_MAX = ((((int64_t)INT_MAX)<<32)|INT_MAX) is the highest
   * supported position.
   *
   * Number of bytes needed for the data is smaller than INT32_MAX
   *
   * qual = NULL, because we do not keep the quality scores through
   * formatting the reads.
   */

  // `qname_nuls` below is the number of '\0' to use to pad the qname
  // so that the cigar has 4-byte alignment.
  const size_t qname_nuls = 4 - l_qname % 4;
  bam_set1_core(bam->core, l_qname, flag, tid, pos, mapq, n_cigar, mtid, mpos,
                isize, l_seq, qname_nuls);

  const size_t data_len = (l_qname + qname_nuls + n_cigar * sizeof(uint32_t) +
                           (l_seq + 1) / 2 + l_seq);

  bam->l_data = data_len;
  if (data_len + l_aux > bam->m_data) {
    const int ret = sam_realloc_bam_data(bam, data_len + l_aux);
    if (ret < 0)
      throw dnmt_error(ret, "Failed to allocate memory for BAM record");
  }
  auto data_iter = bam->data;

  std::copy_n(qname, l_qname, data_iter);
  std::fill_n(data_iter + l_qname, qname_nuls, '\0');
  data_iter += l_qname + qname_nuls;

  // ADS: reinterpret here because we know the cigar is originally an
  // array of uint32_t and has been aligned for efficiency
  std::copy_n(cigar, n_cigar, reinterpret_cast<uint32_t *>(data_iter));
  data_iter += n_cigar * sizeof(uint32_t);

  // skipping sequece assignment
  data_iter += (l_seq + 1) / 2;

  std::fill(data_iter, data_iter + l_seq, '\xff');

  return static_cast<int>(data_len);
}

// static inline size_t
// bam_get_n_cigar(const bam1_t *b) {
//   return b->core.n_cigar;
// }

static inline uint32_t
to_insertion(const uint32_t x) {
  return (x & ~BAM_CIGAR_MASK) | BAM_CINS;
}

static inline void
fix_internal_softclip(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 3) return;
  // find first non-softclip
  auto c_beg = cigar;
  auto c_end = cigar + n_cigar;

  while (!cigar_eats_ref(*c_beg) && ++c_beg != c_end)
    ;
  if (c_beg == c_end) throw dnmt_error("cigar eats no ref");

  while (!cigar_eats_ref(*(c_end - 1)) && --c_end != c_beg)
    ;
  if (c_beg == c_end) throw dnmt_error("cigar eats no ref");

  for (auto c_itr = c_beg; c_itr != c_end; ++c_itr)
    if (bam_cigar_op(*c_itr) == BAM_CSOFT_CLIP) *c_itr = to_insertion(*c_itr);
}

static inline uint32_t
to_softclip(const uint32_t x) {
  return (x & ~BAM_CIGAR_MASK) | BAM_CSOFT_CLIP;
}

static inline void
fix_external_insertion(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 2) return;

  auto c_itr = cigar;
  const auto c_end = c_itr + n_cigar;

  for (; !cigar_eats_ref(*c_itr) && c_itr != c_end; ++c_itr)
    *c_itr = to_softclip(*c_itr);

  if (c_itr == c_end) throw dnmt_error("cigar eats no ref");

  c_itr = cigar + n_cigar - 1;
  for (; !cigar_eats_ref(*c_itr) && c_itr != cigar; --c_itr)
    *c_itr = to_softclip(*c_itr);
}

static inline size_t
merge_cigar_ops(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 2) return n_cigar;
  auto c_itr1 = cigar;
  auto c_end = c_itr1 + n_cigar;
  auto c_itr2 = c_itr1 + 1;
  auto op1 = bam_cigar_op(*c_itr1);
  while (c_itr2 != c_end) {
    auto op2 = bam_cigar_op(*c_itr2);
    if (op1 == op2) {
      *c_itr1 =
        bam_cigar_gen(bam_cigar_oplen(*c_itr1) + bam_cigar_oplen(*c_itr2), op1);
    }
    else {
      *(++c_itr1) = *c_itr2;
      op1 = op2;
    }
    ++c_itr2;
  }
  // another increment to move past final "active" element for c_itr1
  ++c_itr1;
  return std::distance(cigar, c_itr1);
}

static inline size_t
correct_cigar(bam1_t *b) {
  /* This function will change external insertions into soft clip
     operations. Not sure why those would be present. It will also
     change internal soft-clip operations into insertions. This could
     be needed if soft-clipped ends of reads were moved to the middle
     of a merged fragment. Finally, it will collapse adjacent
     identical operations. None of this impacts the seq/qual/aux which
     get moved as a block */

  uint32_t *cigar = bam_get_cigar(b);
  size_t n_cigar = b->core.n_cigar;
  fix_external_insertion(n_cigar, cigar);
  fix_internal_softclip(n_cigar, cigar);

  // merge identical adjacent cigar ops and get new number of ops
  n_cigar = merge_cigar_ops(n_cigar, cigar);
  // difference in bytes to shift the internal data
  const size_t delta = (b->core.n_cigar - n_cigar) * sizeof(uint32_t);
  if (delta > 0) {  // if there is a difference; do the shift
    const auto data_end =
      b->data + b->l_data;  // bam_get_aux(b) + bam_get_l_aux(b);
    std::copy(bam_get_seq(b), data_end, bam_get_seq(b) - delta);
    b->core.n_cigar = n_cigar;  // and update number of cigar ops
  }
  return delta;
}

size_t
correct_cigar(bam_rec &b) {
  return (b.b) ? correct_cigar(b.b) : 0ul;
}

static inline hts_pos_t
get_rlen(const bam1_t *b) {  // less tedious
  return bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
}

static inline size_t
get_l_qseq(const bam1_t *b) {
  return b->core.l_qseq;
}

// static inline void
// complement_seq(char *first, char *last) {
//   for (; first != last; ++first) {
//     assert(valid_base(*first));
//     *first = complement(*first);
//   }
// }

// static inline void
// reverse(char *a, char *b) {
//   char *p1, *p2;
//   for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
//     *p1 ^= *p2;
//     *p2 ^= *p1;
//     *p1 ^= *p2;
//     assert(valid_base(*p1) && valid_base(*p2));
//   }
// }

// return value is the number of cigar ops that are fully consumed in
// order to read n_ref, while "partial_oplen" is the number of bases
// that could be taken from the next operation, which might be merged
// with the other read.
static inline uint32_t
get_full_and_partial_ops(const uint32_t *cig_in, const uint32_t in_ops,
                         const uint32_t n_ref_full, uint32_t *partial_oplen) {
  // assume: n_ops <= size(cig_in) <= size(cig_out)
  size_t rlen = 0;
  uint32_t i = 0;
  for (i = 0; i < in_ops; ++i) {
    if (cigar_eats_ref(cig_in[i])) {
      if (rlen + bam_cigar_oplen(cig_in[i]) > n_ref_full) break;
      rlen += bam_cigar_oplen(cig_in[i]);
    }
  }
  *partial_oplen = n_ref_full - rlen;
  return i;
}

/* This table converts 2 bases packed in a byte to their reverse
 * complement. The input is therefore a unit8_t representing 2 bases.
 * It is assumed that the input uint8_t value is of form "xx" or "x-",
 * where 'x' a 4-bit number representing either A, C, G, T, or N and
 * '-' is 0000.  For example, the ouptut for "AG" is "CT". The format
 * "x-" is often used at the end of an odd-length sequence.  The
 * output of "A-" is "-T", and the output of "C-" is "-G", and so
 * forth. The user must handle this case separately.
 */
const uint8_t byte_revcom_table[] = {
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   8,  136,
  72, 0, 40, 0,   0,  0, 24, 0,   0,  0, 0,  0,   0,  248, 4,  132, 68, 0,
  36, 0, 0,  0,   20, 0, 0,  0,   0,  0, 0,  244, 0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 2,  130, 66, 0,   34, 0,   0,  0,
  18, 0, 0,  0,   0,  0, 0,  242, 0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 1,  129, 65, 0, 33, 0,   0,  0, 17, 0,   0,  0,   0,  0,   0,  241,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 0,  0,   0,  0, 0,  0,   0,  0,   0,  0,   0,  0,
  0,  0, 0,  0,   0,  0, 15, 143, 79, 0, 47, 0,   0,  0,   31, 0,   0,  0,
  0,  0, 0,  255};

static inline void
revcom_byte_then_reverse(unsigned char *a, unsigned char *b) {
  unsigned char *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 = byte_revcom_table[*p1];
    *p2 = byte_revcom_table[*p2];
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
  }
  if (p1 == p2) *p1 = byte_revcom_table[*p1];
}

static inline void
revcomp_seq_by_byte(bam1_t *aln) {
  const size_t l_qseq = get_l_qseq(aln);
  auto seq = bam_get_seq(aln);
  const size_t num_bytes = ceil(l_qseq / 2.0);
  auto seq_end = seq + num_bytes;
  revcom_byte_then_reverse(seq, seq_end);
  if (l_qseq % 2 == 1) {  // for odd-length sequences
    for (size_t i = 0; i < num_bytes - 1; i++) {
      // swap 4-bit chunks within consecutive bytes like this:
      // (----aaaa bbbbcccc dddd....) => (aaaabbbb ccccdddd ....)
      seq[i] = (seq[i] << 4) | (seq[i + 1] >> 4);
    }
    seq[num_bytes - 1] <<= 4;
  }
}

// places seq of b at the end of seq of c
// assumes 0 < c_seq_len - b_seq_len <= a_seq_len
// also assumes that c_seq_len has been figured out
// Also assumes the number of bytes allocated to sequence potion of c->data
// has been set to ceil((a_used_len + b_seq_len) / 2.0) where
// a_used_len = c_seq_len - b_seq_len
static inline void
merge_by_byte(const bam1_t *a, const bam1_t *b, bam1_t *c) {
  // ADS: (todo) need some functions for int_ceil and is_odd
  const size_t b_seq_len = get_l_qseq(b);
  const size_t c_seq_len = get_l_qseq(c);
  const size_t a_used_len = c_seq_len - b_seq_len;

  const bool is_a_odd = a_used_len % 2 == 1;
  const bool is_b_odd = b_seq_len % 2 == 1;
  const bool is_c_odd = c_seq_len % 2 == 1;

  const size_t a_num_bytes = ceil(a_used_len / 2.0);
  const size_t b_num_bytes = ceil(b_seq_len / 2.0);

  const size_t b_offset = is_a_odd && is_b_odd;

  const auto a_seq = bam_get_seq(a);
  const auto b_seq = bam_get_seq(b);
  auto c_seq = bam_get_seq(c);

  std::copy_n(a_seq, a_num_bytes, c_seq);
  if (is_a_odd) {
    // c_seq looks like either [ aa aa aa aa ]
    //                      or [ aa aa aa a- ]
    c_seq[a_num_bytes - 1] &= 0xf0;
    c_seq[a_num_bytes - 1] |=
      is_b_odd ? byte_revcom_table[b_seq[b_num_bytes - 1]]
               : byte_revcom_table[b_seq[b_num_bytes - 1]] >> 4;
  }
  if (is_c_odd) {
    // c_seq looks like either [ aa aa aa aa ]
    //                      or [ aa aa aa ab ]
    for (size_t i = 0; i < b_num_bytes - 1; i++) {
      c_seq[a_num_bytes + i] =
        (byte_revcom_table[b_seq[b_num_bytes - i - 1]] << 4) |
        (byte_revcom_table[b_seq[b_num_bytes - i - 2]] >> 4);
    }
    c_seq[a_num_bytes + b_num_bytes - 1] = byte_revcom_table[b_seq[0]] << 4;
    // Here, c_seq is either [ aa aa aa aa bb bb bb b- ] (a even; b odd)
    //                    or [ aa aa aa ab bb bb bb b- ] (a odd; b odd)
  }
  else {
    for (size_t i = 0; i < b_num_bytes - b_offset; i++) {
      c_seq[a_num_bytes + i] =
        byte_revcom_table[b_seq[b_num_bytes - i - 1 - b_offset]];
    }
    // Here, c_seq is either [ aa aa aa aa bb bb bb bb ] (a even and b even)
    //                    or [ aa aa aa ab bb bb bb    ] (a odd and b odd)
  }
}

static inline void
flip_conversion(bam1_t *aln) {
  aln->core.flag ^= BAM_FREVERSE;  // ADS: flip the "reverse" bit

  revcomp_seq_by_byte(aln);

  // ADS: don't like *(cv + 1) below, but no HTSlib function for it?
  auto cv = bam_aux_get(aln, "CV");
  if (!cv) throw dnmt_error("bam_aux_get failed for CV");
  *(cv + 1) = 'T';
}

void
flip_conversion(bam_rec &aln) {
  flip_conversion(aln.b);
}

// static inline bool
// are_mates(const bam1_t *one, const bam1_t *two) {
//   return one->core.mtid == two->core.tid && one->core.mpos == two->core.pos &&
//          (one->core.flag & BAM_FREVERSE) != (one->core.flag & BAM_FREVERSE);
//   // below is a consistency check and should not be necessary
//   /* &&
//      two->core.mtid == one->core.tid &&
//      two->core.mpos == one->core.pos; */
// }

static inline int
truncate_overlap(const bam1_t *a, const uint32_t overlap, bam1_t *c) {
  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;

  uint32_t part_op = 0;
  const uint32_t c_cur =
    get_full_and_partial_ops(a_cig, a_ops, overlap, &part_op);

  // ADS: hack here because the get_full_and_partial_ops doesn't do
  // exactly what is needed for this.
  const bool use_partial = (c_cur < a->core.n_cigar && part_op > 0);

  const uint32_t c_ops = c_cur + use_partial;
  vector<uint32_t> c_cig(c_ops, 0);

  // ADS: replace this with a std::copy
  auto c_cig_itr = std::copy(a_cig, a_cig + c_cur, begin(c_cig));
  // ADS: warning, if (use_partial==false), then the amount of part_op
  // used below would make no sense.
  if (use_partial)
    *c_cig_itr = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));
  /* after this point the cigar is set and should decide everything */

  const uint32_t c_seq_len = bam_cigar2qlen(c_ops, c_cig.data());
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig.data());

  // flag only needs to worry about strand and single-end stuff
  const uint16_t flag = a->core.flag & (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);

  const size_t a_qname_len = a->core.l_qname - (a->core.l_extranul + 1);
  int ret = bam_set1_wrapper(c, a_qname_len, bam_get_qname(a),
                             flag,  // flags (SR and revcomp info)
                             a->core.tid, a->core.pos, a->core.qual,
                             c_ops,         // merged cigar ops
                             c_cig.data(),  // merged cigar
                             -1,            // (no mate)
                             -1,            // (no mate)
                             isize,         // rlen from new cigar
                             c_seq_len,     // truncated seq length
                             8);            // enough for the 2 tags?
  if (ret < 0) throw dnmt_error(ret, "bam_set1_wrapper");

  const size_t n_bytes_to_copy = (c_seq_len + 1) / 2;  // compression
  std::copy_n(bam_get_seq(a), n_bytes_to_copy, bam_get_seq(c));

  /* add the tags */
  const int64_t nm = bam_aux2i(bam_aux_get(a, "NM"));  // ADS: do better here!
  // "_udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_update_int");

  const uint8_t conversion = bam_aux2A(bam_aux_get(a, "CV"));
  // "_append" for "char" because there is no corresponding update
  ret = bam_aux_append(c, "CV", 'A', 1, &conversion);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_append");

  return ret;
}

int
truncate_overlap(const bam_rec &a, const uint32_t overlap, bam_rec &c) {
  if (c.b == nullptr) c.b = bam_init1();
  return truncate_overlap(a.b, overlap, c.b);
}

int
merge_overlap(const bam1_t *a, const bam1_t *b, const uint32_t head,
              bam1_t *c) {
  assert(head > 0);

  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;

  const uint32_t *b_cig = bam_get_cigar(b);
  const uint32_t b_ops = b->core.n_cigar;

  uint32_t part_op = 0;
  uint32_t c_cur = get_full_and_partial_ops(a_cig, a_ops, head, &part_op);
  // ADS: hack here because the get_full_and_partial_ops doesn't do
  // exactly what is needed for this.
  const bool use_partial = (c_cur < a->core.n_cigar && part_op > 0);

  // check if the middle op would be the same
  const bool merge_mid =
    (use_partial > 0
       ? bam_cigar_op(a_cig[c_cur]) == bam_cigar_op(b_cig[0])
       : bam_cigar_op(a_cig[c_cur - 1]) == bam_cigar_op(b_cig[0]));

  // c_ops: include the prefix of a_cig we need; then add for the
  // partial op; subtract for the identical op in the middle; finally
  // add the rest of b_cig.
  const uint32_t c_ops = c_cur + use_partial - merge_mid + b_ops;
  vector<uint32_t> c_cig(c_ops, 0);
  std::copy(a_cig, a_cig + c_cur, begin(c_cig));

  if (use_partial) {
    c_cig[c_cur] = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));
    c_cur++;  // index of dest for copying b_cig; faciltates corner case
  }
  // Here we get the length of a's sequence part contribution to c's
  // sequence before the possibility of merging the last entry with
  // the first entry in b's cigar. This is done with the cigar, so
  // everything depends on the "use_partial"
  const size_t a_seq_len = bam_cigar2qlen(c_cur, c_cig.data());
  /* ADS: above the return type of bam_cigar2qlen is uint64_t, but
     according to the source as of 05/2023 it cannot become
     negative; no possible error code returned */

  if (merge_mid)  // update the middle op if it's the same
    c_cig[c_cur - 1] = bam_cigar_gen(bam_cigar_oplen(c_cig[c_cur - 1]) +
                                       bam_cigar_oplen(b_cig[0]),
                                     bam_cigar_op(b_cig[0]));
  // copy the cigar from b into c
  std::copy(b_cig + merge_mid, b_cig + b_ops, begin(c_cig) + c_cur);
  /* done with cigar here */

  const uint32_t c_seq_len = a_seq_len + b->core.l_qseq;

  // get the template length
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig.data());

  // flag only needs to worry about strand and single-end stuff
  const uint16_t flag = a->core.flag & (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);

  const size_t a_qname_len = a->core.l_qname - (a->core.l_extranul + 1);
  int ret = bam_set1_wrapper(c, a_qname_len, bam_get_qname(a),
                             flag,  // (no PE; revcomp info)
                             a->core.tid, a->core.pos,
                             a->core.qual,  // mapq from "a" (consider update)
                             c_ops,         // merged cigar ops
                             c_cig.data(),  // merged cigar
                             -1,            // (no mate)
                             -1,            // (no mate)
                             isize,         // updated
                             c_seq_len,     // merged sequence length
                             8);            // enough for 2 tags?
  if (ret < 0) throw dnmt_error(ret, "bam_set1_wrapper in merge_overlap");
  // Merge the sequences by bytes
  merge_by_byte(a, b, c);

  // add the tag for mismatches
  const int64_t nm =
    (bam_aux2i(bam_aux_get(a, "NM")) + bam_aux2i(bam_aux_get(b, "NM")));
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_update_int in merge_overlap");

  // add the tag for conversion
  const uint8_t cv = bam_aux2A(bam_aux_get(a, "CV"));
  ret = bam_aux_append(c, "CV", 'A', 1, &cv);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_append in merge_overlap");

  return ret;
}

int
merge_overlap(const bam_rec &a, const bam_rec &b, const uint32_t head,
              bam_rec &c) {
  if (c.b == nullptr) c.b = bam_init1();
  return merge_overlap(a.b, b.b, head, c.b);
}

static inline int
merge_non_overlap(const bam1_t *a, const bam1_t *b, const uint32_t spacer,
                  bam1_t *c) {
  /* make the cigar string */
  // collect info about the cigar strings
  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;
  const uint32_t *b_cig = bam_get_cigar(b);
  const uint32_t b_ops = b->core.n_cigar;

  // allocate the new cigar string
  const uint32_t c_ops = a_ops + b_ops + 1;
  vector<uint32_t> c_cig(c_ops, 0);

  // concatenate the new cigar strings with a "skip" in the middle
  auto c_cig_itr = std::copy(a_cig, a_cig + a_ops, begin(c_cig));
  *c_cig_itr++ = bam_cigar_gen(spacer, BAM_CREF_SKIP);
  std::copy(b_cig, b_cig + b_ops, c_cig_itr);
  /* done with cigars */

  const size_t a_seq_len = get_l_qseq(a);
  const size_t b_seq_len = get_l_qseq(b);
  const size_t c_seq_len = a_seq_len + b_seq_len;

  // get the template length from the cigar
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig.data());

  // flag: only need to keep strand and single-end info
  const uint16_t flag = a->core.flag & (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);

  const size_t a_qname_len = a->core.l_qname - (a->core.l_extranul + 1);
  int ret = bam_set1_wrapper(c, a_qname_len, bam_get_qname(a),
                             flag,  // flags (no PE; revcomp info)
                             a->core.tid, a->core.pos,
                             a->core.qual,  // mapq from a (consider update)
                             c_ops,         // merged cigar ops
                             c_cig.data(),  // merged cigar
                             -1,            // (no mate)
                             -1,            // (no mate)
                             isize,  // TLEN (relative to reference; SAM docs)
                             c_seq_len,  // merged sequence length
                             8);         // enough for 2 tags of 1 byte value?
  if (ret < 0) throw dnmt_error(ret, "bam_set1 in merge_non_overlap");

  merge_by_byte(a, b, c);

  /* add the tags */
  const int64_t nm =
    (bam_aux2i(bam_aux_get(a, "NM")) + bam_aux2i(bam_aux_get(b, "NM")));
  // "udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) throw dnmt_error(ret, "merge_non_overlap:bam_aux_update_int");

  const uint8_t cv = bam_aux2A(bam_aux_get(a, "CV"));
  // "append" for "char" because there is no corresponding update
  ret = bam_aux_append(c, "CV", 'A', 1, &cv);
  if (ret < 0) throw dnmt_error(ret, "merge_non_overlap:bam_aux_append");

  return ret;
}

int
merge_non_overlap(const bam_rec &a, const bam_rec &b, const uint32_t spacer,
                  bam_rec &c) {
  if (c.b == nullptr) c.b = bam_init1();
  return merge_non_overlap(a.b, b.b, spacer, c.b);
}

static inline int
keep_better_end(const bam1_t *a, const bam1_t *b, bam1_t *c) {
  const auto a_rl = get_rlen(a);
  const auto b_rl = get_rlen(b);
  const auto c_rl = std::max(a_rl, b_rl);
  c = bam_copy1(c, a_rl >= b_rl ? a : b);
  c->core.mtid = -1;
  c->core.mpos = -1;
  c->core.isize = c_rl;
  c->core.flag &= (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);
  return 0;
}

int
keep_better_end(const bam_rec &a, const bam_rec &b, bam_rec &c) {
  if (c.b == nullptr) c.b = bam_init1();
  return keep_better_end(a.b, b.b, c.b);
}

// ADS: will move to using this function once it is written
static inline void
standardize_format(const string &input_format, bam1_t *aln) {
  int err_code = 0;

  if (input_format == "abismal" || input_format == "walt") return;

  if (input_format == "bsmap") {
    // A/T rich: get the ZS tag value
    const auto zs_tag = bam_aux_get(aln, "ZS");
    if (!zs_tag) throw dnmt_error("bam_aux_get for ZS (invalid bsmap)");
    // ADS: test for errors on the line below
    const auto zs_tag_value = string(bam_aux2Z(zs_tag));
    if (zs_tag_value.empty()) throw dnmt_error("empty ZS tag in bsmap format");
    if (zs_tag_value[0] != '-' && zs_tag_value[0] != '+')
      throw dnmt_error("invalid ZS tag in bsmap format");
    const uint8_t cv = zs_tag_value[1] == '-' ? 'A' : 'T';
    // get the "mismatches" tag
    const auto nm_tag = bam_aux_get(aln, "NM");
    if (!nm_tag) throw dnmt_error("invalid NM tag in bsmap format");
    const int64_t nm = bam_aux2i(nm_tag);

    // ADS: this should delete the aux data by truncating the used
    // range within the bam1_t while avoiding resizing memory
    aln->l_data = bam_get_aux(aln) - aln->data;

    /* add the tags we want */
    // "udpate" for "int" because it determines the right size; even
    // though we just deleted all tags, it will add it back here
    err_code = bam_aux_update_int(aln, "NM", nm);
    if (err_code < 0)
      throw dnmt_error(err_code, "error setting NM in bsmap format");

    // "append" for "char" because there is no corresponding update
    err_code = bam_aux_append(aln, "CV", 'A', 1, &cv);
    if (err_code < 0)
      throw dnmt_error(err_code, "error setting conversion in bsmap format");

    // reverse complement if needed
    if (bam_is_rev(aln)) revcomp_seq_by_byte(aln);
  }
  else if (input_format == "bismark") {
    // ADS: Previously we modified the read names at the first
    // underscore. Even if the names are still that way, it should no
    // longer be needed since we compare names up to a learned suffix.

    // A/T rich; get the XR tag value
    auto xr_tag = bam_aux_get(aln, "XR");
    if (!xr_tag) throw dnmt_error("bam_aux_get for XR (invalid bismark)");
    const uint8_t cv = string(bam_aux2Z(xr_tag)) == "GA" ? 'A' : 'T';
    // get the "mismatches" tag
    auto nm_tag = bam_aux_get(aln, "NM");
    if (!nm_tag) throw dnmt_error("bam_aux_get for NM (invalid bismark)");
    const int64_t nm = bam_aux2i(nm_tag);

    aln->l_data = bam_get_aux(aln) - aln->data;  // del aux (no data resize)

    /* add the tags we want */
    // "udpate" for "int" because it determines the right size; even
    // though we just deleted all tags, it will add it back here.
    err_code = bam_aux_update_int(aln, "NM", nm);
    if (err_code < 0) throw dnmt_error(err_code, "bam_aux_update_int");
    // "append" for "char" because there is no corresponding update
    err_code = bam_aux_append(aln, "CV", 'A', 1, &cv);
    if (err_code < 0) throw dnmt_error(err_code, "bam_aux_append");

    if (bam_is_rev(aln))
      revcomp_seq_by_byte(aln);  // reverse complement if needed
  }
  // ADS: the condition below should be checked much earlier, ideally
  // before the output file is created
  else throw runtime_error("incorrect format specified: " + input_format);

  // Be sure this doesn't depend on mapper! Removes the "qual" part of
  // the data in a bam1_t struct but does not change its uncompressed
  // size.
  const auto qs = bam_get_qual(aln);
  std::fill(qs, qs + aln->core.l_qseq, '\xff');  // overwrites qseq
}

void
standardize_format(const string &input_format, bam_rec &aln) {
  standardize_format(input_format, aln.b);
}

// used in methstates
void
apply_cigar(const bam_rec &aln, string &to_inflate,
            const char inflation_symbol) {
  string inflated_seq;
  stringstream ss_cigar;
  size_t i = 0;
  auto to_inflate_beg = std::begin(to_inflate);

  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const auto op = bam_cigar_op(*c_itr);
    const auto n = bam_cigar_oplen(*c_itr);
    ss_cigar << n << op;

    if (cigar_eats_ref(op) && cigar_eats_query(op)) {
      inflated_seq.append(to_inflate_beg + i, to_inflate_beg + i + n);
      i += n;
    }
    else if (cigar_eats_query(op)) {
      // no addition of symbols to query
      i += n;
    }
    else if (cigar_eats_ref(op)) {
      inflated_seq.append(n, inflation_symbol);
      // no increment of index within query
    }
  }

  // sum of total M/I/S/=/X/N operations must equal length of seq
  const size_t orig_len = to_inflate.length();
  if (i != orig_len)
    throw runtime_error(
      "inconsistent number of qseq ops in cigar: " + to_inflate + " " +
      ss_cigar.str() + " " + to_string(i) + " " + to_string(orig_len));
  to_inflate.swap(inflated_seq);
}

void
get_seq_str(const bam_rec &aln, string &seq_str) {
  size_t qlen = static_cast<size_t>(get_l_qseq(aln));
  seq_str.resize(qlen);
  auto seq = bam_get_seq(aln);
  for (size_t i = 0; i < qlen; i++) {
    seq_str[i] = seq_nt16_str[bam_seqi(seq, i)];
  }
}

string
to_string(const bam_header &hdr, const bam_rec &aln) {
  kstring_t ks = {0, 0, nullptr};
  int ret = sam_format1(hdr.h, aln.b, &ks);
  if (ret < 0) { runtime_error("Can't format record: " + to_string(hdr, aln)); }
  if (ks.s != nullptr) free(ks.s);
  return string(ks.s);
}
