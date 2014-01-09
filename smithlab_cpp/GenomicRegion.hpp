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

#ifndef GENOMIC_REGION_HPP
#define GENOMIC_REGION_HPP

#include "smithlab_utils.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <tr1/unordered_map>
#include <limits>

typedef unsigned chrom_id_type;

class GenomicRegion;

class SimpleGenomicRegion {
public:
  // Cannonical
  SimpleGenomicRegion() : chrom(assign_chrom("(null)")), start(0), end(0) {}
  void swap(SimpleGenomicRegion &rhs) {
    std::swap(chrom, rhs.chrom);
    std::swap(start, rhs.start);
    std::swap(end, rhs.end);
  }
  SimpleGenomicRegion& operator=(const SimpleGenomicRegion& rhs) {
    SimpleGenomicRegion tmp(rhs);
    swap(tmp);
    return *this;
  }
  SimpleGenomicRegion(const SimpleGenomicRegion &rhs) :
    chrom(rhs.chrom), start(rhs.start), end(rhs.end) {}

  // Other constructors
  SimpleGenomicRegion(std::string c, size_t sta, size_t e) :
    chrom(assign_chrom(c)), start(sta), end(e) {}
  SimpleGenomicRegion(const GenomicRegion &rhs);
  explicit SimpleGenomicRegion(std::string string_representation);
  SimpleGenomicRegion(const char *string_representation, const size_t len);
  std::string tostring() const;
  
  // accessors
  std::string get_chrom() const {return retrieve_chrom(chrom);}
  size_t get_start() const {return start;}
  size_t get_end() const {return end;}
  size_t get_width() const {return (end > start) ? end - start : 0;}
  
  // mutators
  void set_chrom(const std::string& new_chrom) {chrom = assign_chrom(new_chrom);}
  void set_start(size_t new_start) {start = new_start;}
  void set_end(size_t new_end) {end = new_end;}

  // comparison functions
  bool contains(const SimpleGenomicRegion& other) const;
  bool overlaps(const SimpleGenomicRegion& other) const;
  size_t distance(const SimpleGenomicRegion& other) const;
  bool operator<(const SimpleGenomicRegion& rhs) const;
  bool less1(const SimpleGenomicRegion& rhs) const;
  bool operator<=(const SimpleGenomicRegion& rhs) const;
  bool operator==(const SimpleGenomicRegion& rhs) const;
  bool operator!=(const SimpleGenomicRegion& rhs) const;
  
  bool same_chrom(const SimpleGenomicRegion &other) const {
    return chrom == other.chrom;
  }
  
  friend void
  separate_chromosomes(const std::vector<SimpleGenomicRegion>& regions,
		       std::vector<std::vector<SimpleGenomicRegion> >& 
		       separated_by_chrom);
private:

  static chrom_id_type assign_chrom(const std::string &c);
  static std::string retrieve_chrom(chrom_id_type i);
  
  static std::tr1::unordered_map<std::string, chrom_id_type> fw_table_in;
  static std::tr1::unordered_map<chrom_id_type, std::string> fw_table_out;
  
  // std::string chrom;
  chrom_id_type chrom;
  size_t start;
  size_t end;
};

std::ostream& 
operator<<(std::ostream& the_stream, const SimpleGenomicRegion& region);

std::istream& 
operator>>(std::istream& the_stream, SimpleGenomicRegion& region);

class GenomicRegion {
public:
  GenomicRegion() : chrom(assign_chrom("(NULL)")), 
		    name("X"), start(0), end(0), score(0), strand('+') {}
  void swap(GenomicRegion &rhs) {
    std::swap(chrom, rhs.chrom);
    std::swap(name, rhs.name);
    std::swap(start, rhs.start);
    std::swap(end, rhs.end);
    std::swap(score, rhs.score);
    std::swap(strand, rhs.strand);
  }
  GenomicRegion(const GenomicRegion &other) :
    chrom(other.chrom), name(other.name), start(other.start), end(other.end),
    score(other.score), strand(other.strand) {}
  GenomicRegion& operator=(const GenomicRegion& rhs) {
    GenomicRegion tmp(rhs);
    swap(tmp);
    return *this;
  }
  
  // Other constructors
  GenomicRegion(std::string c, size_t sta, size_t e, 
		std::string n, float sc, char str) :
    chrom(assign_chrom(c)), name(n), start(sta), end(e), score(sc), strand(str) {}
  GenomicRegion(std::string c, size_t sta, size_t e) :
    chrom(assign_chrom(c)), name(std::string("X")), 
    start(sta), end(e), score(0.0), strand('+') {}
  explicit GenomicRegion(std::string string_representation);
  GenomicRegion(const char *s, const size_t len);
  GenomicRegion(const SimpleGenomicRegion &other) :
    chrom(assign_chrom(other.get_chrom())), name("(NULL)"),
    start(other.get_start()), end(other.get_end()), score(0), strand('+') {}
  std::string tostring() const;
  
  // accessors
  std::string get_chrom() const {return retrieve_chrom(chrom);}
  size_t get_start() const {return start;}
  size_t get_end() const {return end;}
  size_t get_width() const {return (end > start) ? end - start : 0;}
  std::string get_name() const {return name;}
  float get_score() const {return score;}
  char get_strand() const {return strand;}
  bool pos_strand() const {return (strand == '+');}
  bool neg_strand() const {return (strand == '-');}

  // mutators
  void set_chrom(const std::string& new_chrom) {chrom = assign_chrom(new_chrom);}
  void set_start(size_t new_start) {start = new_start;}
  void set_end(size_t new_end) {end = new_end;}
  void set_name(const std::string& n) {name = n;}
  void set_score(float s) {score = s;}
  void set_strand(char s) {strand = s;}
  
  // comparison functions
  bool contains(const GenomicRegion& other) const;
  bool overlaps(const GenomicRegion& other) const;
  size_t distance(const GenomicRegion& other) const;
  bool operator<(const GenomicRegion& rhs) const;
  bool less1(const GenomicRegion& rhs) const;
  bool operator<=(const GenomicRegion& rhs) const;
  bool operator!=(const GenomicRegion& rhs) const;
  bool operator==(const GenomicRegion& rhs) const;

  bool same_chrom(const GenomicRegion &other) const {
    return chrom == other.chrom;
  }

  friend void
  separate_chromosomes(const std::vector<GenomicRegion>& regions,
		       std::vector<std::vector<GenomicRegion> >& 
		       separated_by_chrom);
  
private:
  
  static chrom_id_type assign_chrom(const std::string &c);
  static std::string retrieve_chrom(chrom_id_type i);
  
  static std::tr1::unordered_map<std::string, chrom_id_type> fw_table_in;
  static std::tr1::unordered_map<chrom_id_type, std::string> fw_table_out;
  
  // std::string chrom;
  chrom_id_type chrom;
  std::string name;
  size_t start;
  size_t end;
  float score;
  char strand;
};


template <class T> 
bool 
score_less(const T &a, const T &b) {
  return a.get_score() < b.get_score();
}
template <class T> 
bool 
score_greater(const T &a, const T &b) {
  return a.get_score() > b.get_score();
}


class GenomicRegionException : public SMITHLABException {
public:
  GenomicRegionException(std::string s = std::string()) : SMITHLABException(s) {}
};

std::ostream&
operator<<(std::ostream& the_stream, const GenomicRegion& region);

std::istream&
operator>>(std::istream& the_stream, GenomicRegion& region);


template <class T, class U> 
void
sync_chroms(const std::vector<std::vector<T> > &stable, 
	    std::vector<std::vector<U> > &to_sync) {
  std::tr1::unordered_map<std::string, size_t> chrom_index;
  for (size_t i = 0; i < stable.size(); ++i)
    if (!stable[i].empty())
      chrom_index[stable[i].front().get_chrom()] = i;
  std::vector<std::vector<U> > syncd(stable.size());
  for (size_t i = 0; i < to_sync.size(); ++i) {
    if (!to_sync[i].empty()) {
      std::tr1::unordered_map<std::string, size_t>::const_iterator j = 
	chrom_index.find(to_sync[i].front().get_chrom());
      if (j != chrom_index.end())
	to_sync[i].swap(syncd[j->second]);
    }
  }
  syncd.swap(to_sync);
}


template <class T, class U> void
separate_regions(const std::vector<T> &big_regions,
		 const std::vector<U> &regions, 
		 std::vector<std::vector<U> > &sep_regions) {
  size_t rr_id = 0;
  const size_t n_regions = regions.size();
  const size_t n_big_regions = big_regions.size();
  sep_regions.resize(n_big_regions);
  for (size_t i = 0; i < n_big_regions; ++i) {
    const std::string current_chrom(big_regions[i].get_chrom());
    const size_t current_start = big_regions[i].get_start();
    const size_t current_end = big_regions[i].get_end();
    while (rr_id < n_regions &&
	   (regions[rr_id].get_chrom() < current_chrom ||
	    (regions[rr_id].get_chrom() == current_chrom &&
	     regions[rr_id].get_start() < current_start)))
      ++rr_id;
    while (rr_id < n_regions &&
	   (regions[rr_id].get_chrom() == current_chrom &&
	    regions[rr_id].get_start() < current_end)) {
      sep_regions[i].push_back(regions[rr_id]);
      ++rr_id;
    }
  }
}


template <class T> bool 
check_sorted(const std::vector<T> &regions, bool require_unique = false) {
  if (require_unique) {
    for (size_t i = 1; i < regions.size(); ++i)
      if (regions[i] <= regions[i - 1])
	return false;
  }
  else 
    for (size_t i = 1; i < regions.size(); ++i)
      if (regions[i] < regions[i - 1])
	return false;
  return true;
}


template <class T> 
typename std::vector<T>::const_iterator
find_closest(const std::vector<T>& regions, const T& region) {
  typename std::vector<T>::const_iterator closest =
    lower_bound(regions.begin(), regions.end(), region);
  if (closest == regions.begin()) return closest;
  if (closest == regions.end()) return (closest - 1);
  return (region.distance(*closest) < region.distance(*(closest - 1))) ?
    closest : (closest - 1);
}


template <class T> 
typename std::vector<T>::iterator
find_closest(std::vector<T>& regions, const T& region) {
  typename std::vector<T>::iterator closest =
    lower_bound(regions.begin(), regions.end(), region);
  if (closest == regions.begin()) return closest;
  if (closest == regions.end()) return (closest - 1);
  return (region.distance(*closest) < region.distance(*(closest - 1))) ?
    closest : (closest - 1);
}


template <class T> void
collapse(std::vector<T>& regions) {
  typename std::vector<T>::iterator i, good = regions.begin();
  for (i = regions.begin() + 1; i != regions.end(); ++i)
    if (i->overlaps(*good)) {
      good->set_start(std::min(i->get_start(), good->get_start()));
      good->set_end(std::max(i->get_end(), good->get_end()));
    }
    else *(++good) = *i;
  regions.erase(++good, regions.end());
}


template <class T> T
genomic_region_intersection(const T& a, const T& b) {
  if (!a.overlaps(b)) return T(a.get_chrom(), 0, 0);
  else if (a.contains(b)) return b;
  else if (b.contains(a)) return a;
  else if (a < b) return T(a.get_chrom(), b.get_start(), a.get_end());
  else return T(a.get_chrom(), a.get_start(), b.get_end());
}


template <class T> 
void
genomic_region_intersection(const std::vector<T>& regions_a, 
			    const std::vector<T>& regions_b,
			    std::vector<T>& regions_c) {
  typename std::vector<T>::const_iterator a(regions_a.begin());
  typename std::vector<T>::const_iterator a_lim(regions_a.end());
  typename std::vector<T>::const_iterator b(regions_b.begin());
  typename std::vector<T>::const_iterator b_lim(regions_b.end());
  while (a != a_lim && b != b_lim) {
    if (a->overlaps(*b))
      regions_c.push_back(*b);
    if (a == b) {++a; ++b;}
    else if (*a < *b) ++a;
    else ++b; //  if (*b < *a)
  }
}


template <class T> 
void
genomic_region_intersection_by_base(const std::vector<T>& regions_a, 
				    const std::vector<T>& regions_b,
				    std::vector<T>& regions_c) {
  typename std::vector<T>::const_iterator a(regions_a.begin());
  typename std::vector<T>::const_iterator a_lim(regions_a.end());
  typename std::vector<T>::const_iterator b(regions_b.begin());
  typename std::vector<T>::const_iterator b_lim(regions_b.end());
  while (a != a_lim && b != b_lim) {
    if (a->overlaps(*b))
      regions_c.push_back(T(a->get_chrom(), 
			    std::max(a->get_start(), b->get_start()), 
			    std::min(a->get_end(), b->get_end())));
    if (a == b) {++a; ++b;}
    else if (*a < *b) ++a;
    else ++b; //  if (*b < *a)
  }
}

void
ReadBEDFile(std::string filename, std::vector<GenomicRegion> &regions);
void
ReadBEDFile(std::string filename, std::vector<SimpleGenomicRegion> &regions);

template <class T> void
WriteBEDFile(const std::string filename, 
	     const std::vector<std::vector<T> > &regions, 
	     std::string track_name = "") {
  std::ofstream out(filename.c_str());
  if (track_name.length() > 0)
    out << "track name=" << track_name << std::endl;
  for (typename std::vector<std::vector<T> >::const_iterator i = 
	 regions.begin(); i != regions.end(); ++i)
    std::copy(i->begin(), i->end(), std::ostream_iterator<T>(out, "\n"));
  out.close();
}

template <class T> void
WriteBEDFile(const std::string filename,
	     const std::vector<T> &regions, std::string track_name = "") {
  std::ofstream out(filename.c_str());
  if (track_name.length() > 0)
    out << "track name=" << track_name << std::endl;
  std::copy(regions.begin(), regions.end(), std::ostream_iterator<T>(out, "\n"));
  out.close();
}

class BEDFileException : public SMITHLABException {
public:
  BEDFileException(std::string s = std::string()) throw() : SMITHLABException(s) {}
};

void
parse_region_name(std::string region_name,
		  std::string& chrom, size_t &start, size_t &end);

template <class T>
std::string
assemble_region_name(const T &region) {
  return (region.get_chrom() + ":" + smithlab::toa(region.get_start()) + "-" +
          smithlab::toa(region.get_end()));
}

template <class T>
std::string
assemble_region_name(const T &region, const std::string sep) {
  return (region.get_chrom() + sep + smithlab::toa(region.get_start()) + sep +
	  smithlab::toa(region.get_end()));
}

#endif
