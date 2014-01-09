/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
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

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "QualityScore.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <cstring>
#include <cmath>

using std::string;
using std::vector;
using std::ios_base;
using std::cout;

string strip_path(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  return full_path.substr(start);
}

string strip_path_and_suffix(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  size_t end = full_path.find_last_of('.');
  if (end == string::npos)
    end = full_path.length();
  return full_path.substr(start, end - start);
}

void parse_dir_baseanme_suffix(string full_path, string &dirname,
    string &base_name, string &suffix) {
  size_t base_index = full_path.find_last_of("/\\");
  size_t suffix_index = full_path.find_last_of(".");
  if (suffix_index <= base_index)
    suffix_index = full_path.length() - 1;
  dirname = full_path.substr(0, base_index + 1);
  base_name = full_path.substr(base_index + 1, suffix_index - base_index - 1);
  if (suffix_index == full_path.length() - 1)
    suffix = "";
  else
    suffix = full_path.substr(
        suffix_index + 1, full_path.length() - 1 - suffix_index);
}

bool isdir(const char *filename) {
  struct stat buffer;
  stat(filename, &buffer);
  return S_ISDIR(buffer.st_mode);
}

bool is_fastq(const string filename) {
  std::ifstream f(filename.c_str());
  char c = '\0';
  f >> c;
  f.close();
  return (c == '@');
}

////////////////////////////////////////////////////////////////////////
// Stuff dealing with FASTA format sequence files

bool is_valid_filename(const string name, const string& filename_suffix) {
  const string suffix(name.substr(name.find_last_of(".") + 1));
  return (suffix == filename_suffix);
}

string path_join(const string& a, const string& b) {
  if (b.empty() || b[0] == '/')
    throw SMITHLABException("cannot prepend dir to file \"" + b + "\"");
  if (!a.empty() && a[a.length() - 1] == '/')
    return a + b;
  else
    return a + "/" + b;
}

void read_dir(const string& dirname, string filename_suffix,
    vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw SMITHLABException("could not open directory: " + dirname);

  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    if (is_valid_filename(ent->d_name, filename_suffix))
      filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw SMITHLABException("error reading directory: " + dirname);
  if (filenames.empty())
    throw SMITHLABException("no valid files found in: " + dirname);
  closedir(dir);
}

bool is_sequence_line(const char *buffer) {
  return isvalid(buffer[0]);
}

void parse_score_line(const char *buffer, vector<char> &scr) {
  for (const char *i = buffer; *i != '\0'; ++i)
    scr.push_back(*i);
}

inline bool is_fastq_name_line(size_t line_count) {
  return ((line_count & 3ul) == 0ul);
}

inline bool is_fastq_sequence_line(size_t line_count) {
  return ((line_count & 3ul) == 1ul);
}

inline bool is_fastq_score_name_line(size_t line_count) {
  return ((line_count & 3ul) == 2ul);
}

inline bool is_fastq_score_line(size_t line_count) {
  return ((line_count & 3ul) == 3ul);
}

void read_fastq_file(const char *filename, vector<string> &names,
    vector<string> &sequences, vector<vector<double> > &scores) {

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  std::ifstream in(filename);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(filename));

  string s, name;
  vector<char> scr;
  vector<vector<char> > scrs;
  bool first_line = true;
  bool is_sequence_line = false, is_score_line = false;
  size_t line_count = 0;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException(
          "Line in " + name + "\nexceeds max length: "
              + toa(INPUT_BUFFER_SIZE));
    if (in.gcount() == 0)
      break;

    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';

    if (is_fastq_name_line(line_count)) {
      if (buffer[0] != '@')
        throw SMITHLABException("invalid FASTQ name line: " + string(buffer));
      if (first_line == false && s.length() > 0) {
        names.push_back(name);
        sequences.push_back(s);
        scrs.push_back(scr);
      } else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("@ "));
      s = "";
      scr.clear();
      is_sequence_line = true;
    }
    if (is_fastq_sequence_line(line_count)) {
      assert(is_sequence_line);
      s += buffer;
      is_sequence_line = false;
    }
    if (is_fastq_score_name_line(line_count)) {
      if (buffer[0] != '+')
        throw SMITHLABException(
            "invalid FASTQ score name line: " + string(buffer));
      is_score_line = true;
    }
    if (is_fastq_score_line(line_count)) {
      assert(is_score_line);
      parse_score_line(buffer, scr);
      is_score_line = false;
    }
    ++line_count;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
    scrs.push_back(scr);
  }

  using std::ptr_fun;
  using std::not1;
  bool phred_scores = true, solexa_scores = true;
  for (size_t i = 0; i < scrs.size() && phred_scores && solexa_scores; ++i) {
    phred_scores = (phred_scores
        && (find_if(
            scrs[i].begin(), scrs[i].end(), not1(ptr_fun(&valid_phred_score)))
            == scrs[i].end()));
    solexa_scores = (solexa_scores
        && (find_if(
            scrs[i].begin(), scrs[i].end(), not1(ptr_fun(&valid_solexa_score)))
            == scrs[i].end()));
  }

  if (!phred_scores && !solexa_scores)
    throw SMITHLABException("invalid quality scores in FASTQ file");

  for (size_t i = 0; i < scrs.size(); ++i) {
    scores.push_back(vector<double>(scrs[i].size()));
    for (size_t j = 0; j < scrs[i].size(); ++j)
      scores[i][j] =
          (solexa_scores) ?
              quality_character_to_solexa(scrs[i][j] - 5) :
              quality_character_to_phred(scrs[i][j]);
    scrs[i].clear();
  }
}

void read_fastq_file(const char *filename, vector<string> &names,
    vector<string> &sequences, vector<string> &scores) {

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  std::ifstream in(filename);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(filename));

  string s, name, scr;
  bool first_line = true;
  bool is_sequence_line = false, is_score_line = false;
  size_t line_count = 0;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException(
          "Line in " + name + "\nexceeds max length: "
              + toa(INPUT_BUFFER_SIZE));
    if (in.gcount() == 0)
      break;

    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';

    if (is_fastq_name_line(line_count)) {
      if (buffer[0] != '@')
        throw SMITHLABException("invalid FASTQ name line: " + string(buffer));
      if (first_line == false && s.length() > 0) {
        names.push_back(name);
        sequences.push_back(s);
        scores.push_back(scr);
      } else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("@ "));
      is_sequence_line = true;
    }
    if (is_fastq_sequence_line(line_count)) {
      assert(is_sequence_line);
      s = buffer;
      is_sequence_line = false;
    }
    if (is_fastq_score_name_line(line_count)) {
      if (buffer[0] != '+')
        throw SMITHLABException(
            "invalid FASTQ score name line: " + string(buffer));
      is_score_line = true;
    }
    if (is_fastq_score_line(line_count)) {
      assert(is_score_line);
      scr = buffer;
      is_score_line = false;
    }
    ++line_count;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
    scores.push_back(scr);
  }
}

void read_fasta_file(const string filename, vector<string> &names,
    vector<string> &sequences) {

  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) {
    throw SMITHLABException("cannot open input file " + string(filename));
  }

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  string s, name;

  bool first_line = true;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException(
          "Line in " + name + "\nexceeds max length: "
              + toa(INPUT_BUFFER_SIZE));
    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    if (buffer[0] == '>') {
      if (first_line == false && s.length() > 0) {
        names.push_back(name);
        sequences.push_back(s);
      } else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("> "));
      const size_t first_whitespace = name.find_first_of(" \t");
      if (first_whitespace != std::string::npos)
        name = name.substr(0, first_whitespace);
      s = "";
    } else
      s += buffer;
    in.peek();
  } //while
  if (!first_line && s.length() > 0) {

    names.push_back(name);
    sequences.push_back(s);
  }
}

void read_fasta_file(const string filename, const string &target,
    string &sequence) {

  // read the sequence with the given name from a fasta file

  sequence = "";

  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) {
    throw SMITHLABException("cannot open input file " + string(filename));
  }

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  string s, name;

  bool first_line = true;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException(
          "Line in " + name + "\nexceeds max length: "
              + toa(INPUT_BUFFER_SIZE));
    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    if (buffer[0] == '>') {
      if (first_line == false && s.length() > 0 && name == target) {
        std::swap(sequence, s);
        in.close();
        return;
      } else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("> "));
      const size_t first_whitespace = name.find_first_of(" \t");
      if (first_whitespace != std::string::npos)
        name = name.substr(0, first_whitespace);
      s = "";
    } else if (name == target)
      s += buffer;
    in.peek();
  } //while
  if (!first_line && s.length() > 0 && name == target) {
    std::swap(sequence, s);
  }
  in.close();
}

void read_filename_file(const char *filename, vector<string> &filenames) {

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  std::ifstream in(filename);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(filename));
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException(
          "Line in " + string(filename) + "\nexceeds max length: "
              + toa(INPUT_BUFFER_SIZE));
    filenames.push_back(buffer);
    in.peek();
  }
}

size_t get_filesize(string filename) {
  std::ifstream f(filename.c_str());
  if (!f.good()) {
    return 0;
  }
  size_t begin_pos = f.tellg();
  f.seekg(0, ios_base::end);
  size_t end_pos = f.tellg();
  f.close();
  return end_pos - begin_pos;
}

string basename(string filename) {
  const string s(filename.substr(0, filename.find_last_of(".")));
  const size_t final_slash = s.find_last_of("/");
  if (final_slash != string::npos)
    return s.substr(final_slash + 1);
  else
    return s;
}

static size_t adjust_start_pos(const size_t orig_start,
    const string &chrom_name) {
  static const double LINE_WIDTH = 50.0;
  const size_t name_offset = chrom_name.length() + 2; // For the '>' and the '\n';
  const size_t preceding_newlines = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  return orig_start + preceding_newlines + name_offset;
}

static size_t adjust_region_size(const size_t orig_start,
    const string &chrom_name, const size_t orig_size) {
  static const double LINE_WIDTH = 50.0;
  const size_t preceding_newlines_start = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  const size_t preceding_newlines_end = static_cast<size_t>(std::floor(
      (orig_start + orig_size) / LINE_WIDTH));
  return (orig_size + (preceding_newlines_end - preceding_newlines_start));
}

void extract_regions_chrom_fasta(const string &chrom_name,
    const string &filename, const vector<SimpleGenomicRegion> &regions,
    vector<string> &sequences) {

  std::ifstream in(filename.c_str());
  for (vector<SimpleGenomicRegion>::const_iterator i(regions.begin());
      i != regions.end(); ++i) {

    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size = adjust_region_size(
        orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(
        buffer, buffer + region_size,
        std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    std::transform(
        sequences.back().begin(), sequences.back().end(),
        sequences.back().begin(), std::ptr_fun(&toupper));
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}

void extract_regions_chrom_fasta(const string &chrom_name,
    const string &filename, const vector<GenomicRegion> &regions,
    vector<string> &sequences) {

  std::ifstream in(filename.c_str());
  for (vector<GenomicRegion>::const_iterator i(regions.begin());
      i != regions.end(); ++i) {

    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size = adjust_region_size(
        orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(
        buffer, buffer + region_size,
        std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    std::transform(
        sequences.back().begin(), sequences.back().end(),
        sequences.back().begin(), std::ptr_fun(&toupper));
    if (i->neg_strand())
      revcomp_inplace(sequences.back());
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}

void read_dir(const string& dirname, vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw "could not open directory: " + dirname;

  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw "error reading directory: " + dirname;
  if (filenames.empty())
    throw "no valid files found in: " + dirname;
  closedir(dir);
}

void extract_regions_fasta(const string &dirname,
    const vector<GenomicRegion> &regions_in, vector<string> &sequences) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<GenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    // GET THE RIGHT FILE
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
        chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw SMITHLABException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(
        chrom_name, filenames[f_idx->second], regions[i], sequences);
  }
}

void extract_regions_fasta(const string &dirname,
    const vector<SimpleGenomicRegion> &regions_in, vector<string> &sequences) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<SimpleGenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    // GET THE RIGHT FILE
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
        chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw SMITHLABException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(
        chrom_name, filenames[f_idx->second], regions[i], sequences);
  }
}

void read_prb_file(string filename, vector<vector<vector<double> > > &scores) {
  static const size_t INPUT_BUFFER_SIZE = 1000000;
  scores.clear();
  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + filename);
  string s;
  size_t line_number = 0;
  while (!in.eof()) {
    ++line_number;
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException(
          "Line in " + filename + "\nexceeds max length: "
              + toa(INPUT_BUFFER_SIZE));
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';

    vector<string> parts;
    smithlab::split_whitespace(buffer, parts);
    if (parts.size() % smithlab::alphabet_size != 0)
      throw SMITHLABException(
          "Incorrect number of values on line " + toa(line_number) + " in file "
              + filename);
    scores.push_back(vector<vector<double> >());
    for (size_t i = 0; i < parts.size(); i += smithlab::alphabet_size) {
      scores.back().push_back(vector<double>());
      for (size_t j = 0; j < smithlab::alphabet_size; ++j)
        scores.back().back().push_back(atof(parts[i + j].c_str()));
    }
    in.peek();
  }
}
