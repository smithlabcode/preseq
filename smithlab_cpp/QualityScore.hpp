/*
 *    Part of SMITHLAB software
 *
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

#ifndef QUALITY_SCORE_HPP
#define QUALITY_SCORE_HPP

#include <cmath>
#include <algorithm>
#include <string>

////////////////////////////////////////////////////////////////////////
//
// FUNCTIONS FOR MANIPULATING VALUES REALTED TO QUALITY SCORES
//

typedef size_t FASTQScoreType;

static const FASTQScoreType FASTQ_Solexa = 0;
static const FASTQScoreType FASTQ_Phred = 1;

inline bool
FASTQScoreIsSolexa(const FASTQScoreType t) {return t == FASTQ_Solexa;}

inline bool
FASTQScoreIsPhred(const FASTQScoreType t) {return t == FASTQ_Phred;}

static const double neg_ten_over_log_ten = -4.342944819032517501;

//// CONVERT _FROM_ ERROR PROBABILITIES
inline double
error_probability_to_phred(const double r) {
  return neg_ten_over_log_ten*std::log(r);
}
inline double
error_probability_to_solexa(const double r) {
  return neg_ten_over_log_ten*(std::log(r) - std::log(1.0 - r));
}

//// CONVERT _TO_ ERROR PROBABILITIES
inline double
phred_to_error_probability(const double r) {
  const double h = r/neg_ten_over_log_ten;
  return std::exp(h);
}
inline double
solexa_to_error_probability(const double r) {
  const double s = r/neg_ten_over_log_ten;
  return std::exp(s)/(1.0 + std::exp(s));
}

//// CONVERT _FROM_ QUALITY CHARACTERS (I.E. THE CHARACTERS IN FASTQ FILES)
inline char
quality_character_to_phred(const char c) {
  return char(c - 33);
}
inline char
quality_character_to_solexa(const char c) {
  return char(c - 64);
}

//// CONVERT _TO_ QUALITY CHARACTERS (I.E. THE CHARACTERS IN FASTQ FILES)
inline char
phred_to_quality_character(const double h) {
  return char(std::min(60.0, h) + 33);
}
inline char
solexa_to_quality_character(const double s) {
  return char(std::min(40.0, s) + 64);
}

//// CHECK FOR VALIDITY OF A FASTQ SCORE CHARACTER
inline bool
valid_phred_score(const char c) {
  return (c >= 33 && c <= 93);
}
inline bool
valid_solexa_score(const char c) {
  // return (c >= 64 && c <= 104);
  return (c >= 59 && c <= 104); // to allow for old Solexa format
}

inline double
quality_char_to_error_probability(const FASTQScoreType t,
				  const char c) {
  return (t == FASTQ_Solexa) ?
    solexa_to_error_probability(quality_character_to_solexa(c)) :
    phred_to_error_probability(quality_character_to_phred(c));
}

inline double
char2prob_solexa(const char c)
{
    return 1 - solexa_to_error_probability(quality_character_to_solexa(c));
}

inline double
char2prob_phred(const char c)
{
    return 1 - phred_to_error_probability(quality_character_to_phred(c));
}

inline double
char2err_solexa(const char c)
{
    return solexa_to_error_probability(quality_character_to_solexa(c));
}

inline double
char2err_phred(const char c)
{
    return phred_to_error_probability(quality_character_to_phred(c));
}

inline char
prob2char_solexa(const double prob)
{
    return solexa_to_quality_character(
        error_probability_to_solexa(1 - prob));
}

inline char
prob2char_phred(const double prob)
{
    return phred_to_quality_character(
        error_probability_to_phred(1 - prob));
}

inline char
err2char_solexa(const double err)
{
    return solexa_to_quality_character(
        error_probability_to_solexa(err));
}

inline char
err2char_phred(const double err)
{
    return phred_to_quality_character(
        error_probability_to_phred(err));
}

inline double
quality_score_to_error_probability(const FASTQScoreType t,
				   const double s) {
  return (FASTQScoreIsSolexa(t)) ?
    solexa_to_error_probability(s) : phred_to_error_probability(s);
}

FASTQScoreType
fastq_score_type(const std::string filename);
FASTQScoreType
mapped_reads_score_type(const std::string filename);

#endif
