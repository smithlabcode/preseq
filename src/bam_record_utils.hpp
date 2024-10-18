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

#ifndef SRC_BAM_RECORD_UTILS_HPP_
#define SRC_BAM_RECORD_UTILS_HPP_

/* ADS: need to control all the macros from HTSlib pollution. For
   functions maybe:

   $ gcc -dM -E sam.h | grep "define [a-z]" | awk '{print $2}' |\
       grep "[(]" | awk -v FS="(" '{print "#undef",$1}'

   This gives about 65 symbols that need to be deleted. For the others
   I don't know what to do because some of them have "#define _" which
   means they should be system symbols.
*/

#include <bamxx.hpp>

#include <string>

#ifdef bam_is_rev
#undef bam_is_rev
#endif

inline bool
bam_is_rev(const bamxx::bam_rec &b) {
  return (b.b->core.flag & BAM_FREVERSE) != 0;
}

#ifdef bam_is_mrev
#undef bam_is_mrev
#endif

inline bool
bam_is_mrev(const bamxx::bam_rec &b) {
  return (b.b->core.flag & BAM_FMREVERSE) != 0;
}

#ifdef bam_get_qname
#undef bam_get_qname
#endif

inline char *
bam_get_qname(const bamxx::bam_rec &b) {
  return reinterpret_cast<char *>(b.b->data);
}

#ifdef bam_get_cigar
#undef bam_get_cigar
#endif

inline uint32_t *
bam_get_cigar(const bamxx::bam_rec &b) {
  // start of data + bytes for query/read name
  return reinterpret_cast<uint32_t *>(b.b->data + b.b->core.l_qname);
}

#ifdef bam_get_seq
#undef bam_get_seq
#endif

inline uint8_t *
bam_get_seq(const bamxx::bam_rec &b) {
  // start of data + bytes for cigar + bytes for query/read name
  return b.b->data + b.b->core.l_qname + (b.b->core.n_cigar << 2);
}

#ifdef bam_get_qual
#undef bam_get_qual
#endif

inline uint8_t *
bam_get_qual(const bamxx::bam_rec &b) {
  return b.b->data +                     // start of data
         b.b->core.l_qname +             // bytes for query name
         (b.b->core.n_cigar << 2) +      // bytes for cigar
         ((b.b->core.l_qseq + 1) >> 1);  // bytes for packed query/read
}

#ifdef bam_get_aux
#undef bam_get_aux
#endif

inline uint8_t *
bam_get_aux(const bamxx::bam_rec &b) {
  return b.b->data + b.b->core.l_qname + (b.b->core.n_cigar << 2) +
         ((b.b->core.l_qseq + 1) >> 1) + b.b->core.l_qseq;
}

#ifdef bam_get_l_aux
#undef bam_get_l_aux
#endif

inline int
bam_get_l_aux(const bamxx::bam_rec &b) {
  return b.b->l_data - (b.b->core.l_qname + (b.b->core.n_cigar << 2) +
                        ((b.b->core.l_qseq + 1) >> 1) + b.b->core.l_qseq);
}

#ifdef bam_cigar_op
#undef bam_cigar_op
#endif

inline uint32_t
bam_cigar_op(const uint32_t c) {
  return c & BAM_CIGAR_MASK;
}

#ifdef bam_cigar_oplen
#undef bam_cigar_oplen
#endif

inline uint32_t
bam_cigar_oplen(const uint32_t c) {
  return c >> BAM_CIGAR_SHIFT;
}

inline bool
bam_same_orientation(const bamxx::bam_rec &a, const bamxx::bam_rec &b) {
  return ((a.b->core.flag ^ b.b->core.flag) & BAM_FREVERSE) != 0;
}

int
truncate_overlap(const bamxx::bam_rec &a, const uint32_t overlap,
                 bamxx::bam_rec &c);

int
merge_overlap(const bamxx::bam_rec &a, const bamxx::bam_rec &b,
              const uint32_t head, bamxx::bam_rec &c);

int
merge_non_overlap(const bamxx::bam_rec &a, const bamxx::bam_rec &b,
                  const uint32_t spacer, bamxx::bam_rec &c);

int
keep_better_end(const bamxx::bam_rec &a, const bamxx::bam_rec &b,
                bamxx::bam_rec &c);

size_t
correct_cigar(bamxx::bam_rec &b);

void
flip_conversion(bamxx::bam_rec &aln);

inline bool
is_a_rich(const bamxx::bam_rec &b) {
  return bam_aux2A(bam_aux_get(b.b, "CV")) == 'A';
}

void
standardize_format(const std::string &input_format, bamxx::bam_rec &aln);

void
apply_cigar(const bamxx::bam_rec &aln, std::string &to_inflate,
            const char inflation_symbol);

void
get_seq_str(const bamxx::bam_rec &aln, std::string &seq_str);

inline bool
are_mates(const bamxx::bam_rec &one, const bamxx::bam_rec &two) {
  return one.b->core.mtid == two.b->core.tid &&
         one.b->core.mpos == two.b->core.pos && bam_same_orientation(one, two);
  // below is a consistency check and should not be necessary
  /* &&
     two->core.mtid == one->core.tid &&
     two->core.mpos == one->core.pos; */
}

inline int32_t
get_l_qseq(const bamxx::bam_rec &b) {
  return b.b->core.l_qseq;
}

inline size_t
get_n_targets(const bamxx::bam_header &bh) {
  return bh.h->n_targets;
}

inline std::string
get_qname(const bamxx::bam_rec &b) {
  return bam_get_qname(b);
}

inline int32_t
get_tid(const bamxx::bam_rec &b) {
  return b.b->core.tid;
}

inline hts_pos_t
get_pos(const bamxx::bam_rec &b) {
  return b.b->core.pos;
}

inline int32_t
get_mtid(const bamxx::bam_rec &b) {
  return b.b->core.mtid;
}

inline hts_pos_t
get_mpos(const bamxx::bam_rec &b) {
  return b.b->core.mpos;
}

inline uint32_t
get_n_cigar(const bamxx::bam_rec &b) {
  return b.b->core.n_cigar;
}

inline hts_pos_t
get_endpos(const bamxx::bam_rec &b) {
  return bam_endpos(b.b);
}

inline bool
cigar_eats_ref(const uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 2;
}

inline bool
cigar_eats_query(const uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 1;
}

inline bool
cigar_eats_frag(const uint32_t c) {
  return bam_cigar_op(c) == BAM_CREF_SKIP;
}

inline bool
precedes_by_start(const bamxx::bam_rec &a, const bamxx::bam_rec &b) {
  // assumes a.get_tid() <= b.get_tid()
  return get_tid(a) == get_tid(b) && get_pos(a) < get_pos(b);
}

inline bool
precedes_by_end_and_strand(const bamxx::bam_rec &a, const bamxx::bam_rec &b) {
  const auto end_a = bam_endpos(a.b);
  const auto end_b = bam_endpos(b.b);
  return end_a < end_b || (end_a == end_b && bam_is_rev(a) < bam_is_rev(b));
}

inline bool
equivalent_chrom_and_start(const bamxx::bam_rec &a, const bamxx::bam_rec &b) {
  return a.b->core.pos == b.b->core.pos && a.b->core.tid == b.b->core.tid;
}

inline bool
equivalent_end_and_strand(const bamxx::bam_rec &a, const bamxx::bam_rec &b) {
  return bam_endpos(a.b) == bam_endpos(b.b) && bam_is_rev(a) == bam_is_rev(b);
}

template <typename T>
int
bam_aux_update_int(bamxx::bam_rec &b, const char tag[2], T val) {
  return bam_aux_update_int(b.b, tag, val);
}

inline std::string
sam_hdr_tid2name(const bamxx::bam_header &hdr, const int32_t tid) {
  return std::string(sam_hdr_tid2name(hdr.h, tid));
}

inline uint32_t
sam_hdr_tid2len(const bamxx::bam_header &hdr, const int32_t tid) {
  return sam_hdr_tid2len(hdr.h, tid);
}

inline std::string
sam_hdr_tid2name(const bamxx::bam_header &hdr, const bamxx::bam_rec &aln) {
  return std::string(sam_hdr_tid2name(hdr.h, aln.b->core.tid));
}

std::string
to_string(const bamxx::bam_header &hdr, const bamxx::bam_rec &aln);

inline size_t
rlen_from_cigar(const bamxx::bam_rec &aln) {
  return bam_cigar2rlen(get_n_cigar(aln), bam_get_cigar(aln));
}

#endif  // SRC_BAM_RECORD_UTILS_HPP_
