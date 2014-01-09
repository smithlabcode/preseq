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

#ifndef SMITHLAB_OS_HPP
#define SMITHLAB_OS_HPP

#include <string>
#include <vector>

#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>
#include "GenomicRegion.hpp"

bool
isdir(const char *filename);

bool
is_fastq(const std::string filename);

bool
is_valid_filename(const std::string name, 
		  const std::string& filename_suffix);

std::string 
path_join(const std::string& a, const std::string& b);

void 
read_dir(const std::string& dirname, 
	 std::string filename_suffix,
	 std::vector<std::string> &filenames);

void
read_fasta_file(const std::string filename,
		std::vector<std::string> &names, 
		std::vector<std::string> &sequences);

void
read_fasta_file(const std::string filename,
        const std::string &name,
        std::string &sequence);

void
read_fastq_file(const char *filename, 
		std::vector<std::string> &names, 
		std::vector<std::string> &sequences,
		std::vector<std::vector<double> > &scores);

void
read_fastq_file(const char *filename, 
		std::vector<std::string> &names, 
		std::vector<std::string> &sequences,
		std::vector<std::string> &scores);

void
read_prb_file(std::string filename, 
	      std::vector<std::vector<std::vector<double> > > &scores);

void
read_filename_file(const char *filename, 
		   std::vector<std::string> &filenames);

size_t 
get_filesize(std::string filename);

std::string
basename(std::string filename);

void
parse_dir_baseanme_suffix(std::string full_path,
			  std::string &dirname,
			  std::string &base_name,
			  std::string &suffix);

std::string
strip_path(std::string full_path);

std::string
strip_path_and_suffix(std::string full_path);


void
extract_regions_chrom_fasta(const std::string &chrom_name, 
			    const std::string &filename,
			    const std::vector<GenomicRegion> &regions, 
			    std::vector<std::string> &sequences);


void
extract_regions_chrom_fasta(const std::string &chrom_name, 
			    const std::string &filename,
			    const std::vector<SimpleGenomicRegion> &regions, 
			    std::vector<std::string> &sequences);

void
extract_regions_fasta(const std::string &dirname, 
		      const std::vector<SimpleGenomicRegion> &regions_in, 
		      std::vector<std::string> &sequences);

void
extract_regions_fasta(const std::string &dirname, 
		      const std::vector<GenomicRegion> &regions_in, 
		      std::vector<std::string> &sequences);

void 
read_dir(const std::string& dirname, std::vector<std::string> &filenames);

#endif
