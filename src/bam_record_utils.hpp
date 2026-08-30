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

#include "bamxx/bamxx.hpp"

#include <htslib/sam.h>

#include <cstddef>
#include <cstdint>
#include <string>

#ifdef bam_is_rev
#undef bam_is_rev
#endif

inline auto
bam_is_rev(const bamxx::bam_rec &b) -> bool {
  return (b.b->core.flag & BAM_FREVERSE) != 0;
}

#ifdef bam_is_mrev
#undef bam_is_mrev
#endif

inline auto
bam_is_mrev(const bamxx::bam_rec &b) -> bool {
  return (b.b->core.flag & BAM_FMREVERSE) != 0;
}

#ifdef bam_get_qname
#undef bam_get_qname
#endif

inline auto
bam_get_qname(const bamxx::bam_rec &b) -> char * {
  return reinterpret_cast<char *>(b.b->data);
}

#ifdef bam_get_cigar
#undef bam_get_cigar
#endif

inline auto
bam_get_cigar(const bamxx::bam_rec &b) -> uint32_t * {
  // start of data + bytes for query/read name
  return reinterpret_cast<uint32_t *>(b.b->data + b.b->core.l_qname);
}

#ifdef bam_get_seq
#undef bam_get_seq
#endif

inline auto
bam_get_seq(const bamxx::bam_rec &b) -> uint8_t * {
  // start of data + bytes for cigar + bytes for query/read name
  return b.b->data + b.b->core.l_qname + (b.b->core.n_cigar << 2);
}

#ifdef bam_get_qual
#undef bam_get_qual
#endif

inline auto
bam_get_qual(const bamxx::bam_rec &b) -> uint8_t * {
  return b.b->data +                     // start of data
         b.b->core.l_qname +             // bytes for query name
         (b.b->core.n_cigar << 2) +      // bytes for cigar
         ((b.b->core.l_qseq + 1) >> 1);  // bytes for packed query/read
}

#ifdef bam_get_aux
#undef bam_get_aux
#endif

inline auto
bam_get_aux(const bamxx::bam_rec &b) -> uint8_t * {
  return b.b->data + b.b->core.l_qname + (b.b->core.n_cigar << 2) +
         ((b.b->core.l_qseq + 1) >> 1) + b.b->core.l_qseq;
}

#ifdef bam_get_l_aux
#undef bam_get_l_aux
#endif

inline auto
bam_get_l_aux(const bamxx::bam_rec &b) -> int {
  return b.b->l_data - (b.b->core.l_qname + (b.b->core.n_cigar << 2) +
                        ((b.b->core.l_qseq + 1) >> 1) + b.b->core.l_qseq);
}

#ifdef bam_cigar_op
#undef bam_cigar_op
#endif

inline auto
bam_cigar_op(const uint32_t c) -> uint32_t {
  return c & BAM_CIGAR_MASK;
}

#ifdef bam_cigar_oplen
#undef bam_cigar_oplen
#endif

inline auto
bam_cigar_oplen(const uint32_t c) -> uint32_t {
  return c >> BAM_CIGAR_SHIFT;
}

inline auto
bam_same_orientation(const bamxx::bam_rec &a, const bamxx::bam_rec &b) -> bool {
  return ((a.b->core.flag ^ b.b->core.flag) & BAM_FREVERSE) != 0;
}

auto
truncate_overlap(const bamxx::bam_rec &a, const uint32_t overlap,
                 bamxx::bam_rec &c) -> int;

auto
merge_overlap(const bamxx::bam_rec &a, const bamxx::bam_rec &b,
              const uint32_t head, bamxx::bam_rec &c) -> int;

auto
merge_non_overlap(const bamxx::bam_rec &a, const bamxx::bam_rec &b,
                  const uint32_t spacer, bamxx::bam_rec &c) -> int;

auto
keep_better_end(const bamxx::bam_rec &a, const bamxx::bam_rec &b,
                bamxx::bam_rec &c) -> int;

auto
correct_cigar(bamxx::bam_rec &b) -> std::size_t;

void
flip_conversion(bamxx::bam_rec &aln);

inline auto
is_a_rich(const bamxx::bam_rec &b) -> bool {
  return bam_aux2A(bam_aux_get(b.b, "CV")) == 'A';
}

void
standardize_format(const std::string &input_format, bamxx::bam_rec &aln);

void
apply_cigar(const bamxx::bam_rec &aln, std::string &to_inflate,
            const char inflation_symbol);

void
get_seq_str(const bamxx::bam_rec &aln, std::string &seq_str);

inline auto
are_mates(const bamxx::bam_rec &one, const bamxx::bam_rec &two) -> bool {
  return one.b->core.mtid == two.b->core.tid &&
         one.b->core.mpos == two.b->core.pos && bam_same_orientation(one, two);
  // below is a consistency check and should not be necessary
  /* &&
     two->core.mtid == one->core.tid &&
     two->core.mpos == one->core.pos; */
}

inline auto
get_l_qseq(const bamxx::bam_rec &b) -> int32_t {
  return b.b->core.l_qseq;
}

inline auto
get_n_targets(const bamxx::bam_header &bh) -> std::size_t {
  return bh.h->n_targets;
}

inline auto
get_qname(const bamxx::bam_rec &b) -> std::string {
  return bam_get_qname(b);
}

inline auto
get_tid(const bamxx::bam_rec &b) -> int32_t {
  return b.b->core.tid;
}

inline auto
get_pos(const bamxx::bam_rec &b) -> hts_pos_t {
  return b.b->core.pos;
}

inline auto
get_mtid(const bamxx::bam_rec &b) -> int32_t {
  return b.b->core.mtid;
}

inline auto
get_mpos(const bamxx::bam_rec &b) -> hts_pos_t {
  return b.b->core.mpos;
}

inline auto
get_n_cigar(const bamxx::bam_rec &b) -> uint32_t {
  return b.b->core.n_cigar;
}

inline auto
get_endpos(const bamxx::bam_rec &b) -> hts_pos_t {
  return bam_endpos(b.b);
}

inline auto
cigar_eats_ref(const uint32_t c) -> bool {
  return bam_cigar_type(bam_cigar_op(c)) & 2;
}

inline auto
cigar_eats_query(const uint32_t c) -> bool {
  return bam_cigar_type(bam_cigar_op(c)) & 1;
}

inline auto
cigar_eats_frag(const uint32_t c) -> bool {
  return bam_cigar_op(c) == BAM_CREF_SKIP;
}

inline auto
precedes_by_start(const bamxx::bam_rec &a, const bamxx::bam_rec &b) -> bool {
  // assumes a.get_tid() <= b.get_tid()
  return get_tid(a) == get_tid(b) && get_pos(a) < get_pos(b);
}

inline auto
precedes_by_end_and_strand(const bamxx::bam_rec &a,
                           const bamxx::bam_rec &b) -> bool {
  const auto end_a = bam_endpos(a.b);
  const auto end_b = bam_endpos(b.b);
  return end_a < end_b ||
         (end_a == end_b && bam_is_rev(a) == false && bam_is_rev(b) == true);
}

inline auto
equivalent_chrom_and_start(const bamxx::bam_rec &a,
                           const bamxx::bam_rec &b) -> bool {
  return a.b->core.pos == b.b->core.pos && a.b->core.tid == b.b->core.tid;
}

inline auto
equivalent_end_and_strand(const bamxx::bam_rec &a,
                          const bamxx::bam_rec &b) -> bool {
  return bam_endpos(a.b) == bam_endpos(b.b) && bam_is_rev(a) == bam_is_rev(b);
}

template <typename T>
auto
bam_aux_update_int(bamxx::bam_rec &b, const char tag[2], T val) -> int {
  return bam_aux_update_int(b.b, tag, val);
}

inline auto
sam_hdr_tid2name(const bamxx::bam_header &hdr,
                 const int32_t tid) -> std::string {
  return std::string(sam_hdr_tid2name(hdr.h, tid));
}

inline auto
sam_hdr_tid2len(const bamxx::bam_header &hdr, const int32_t tid) -> uint32_t {
  return sam_hdr_tid2len(hdr.h, tid);
}

inline auto
sam_hdr_tid2name(const bamxx::bam_header &hdr,
                 const bamxx::bam_rec &aln) -> std::string {
  return std::string(sam_hdr_tid2name(hdr.h, aln.b->core.tid));
}

auto
to_string(const bamxx::bam_header &hdr,
          const bamxx::bam_rec &aln) -> std::string;

inline auto
rlen_from_cigar(const bamxx::bam_rec &aln) -> std::size_t {
  return bam_cigar2rlen(get_n_cigar(aln), bam_get_cigar(aln));
}

#endif  // SRC_BAM_RECORD_UTILS_HPP_
