#preseq

The preseq package is aimed at predicting the yield of distinct 
reads from a genomic library from an initial sequencing 
experiment. The estimates can then be used to examine the utility 
of further sequencing, optimize the sequencing depth, or to 
screen multiple libraries to avoid low complexity samples.

preseq is fully documented in the manual: [download here](https://github.com/smithlabcode/preseq/raw/master/docs/manual.pdf)

## Readme Table of Contents
* [Updates from Previous Releases](#updates-from-previous-releases)
* [Contact Information](#contact-information)
* [System Requirements](#system-requirements)
* [Installation](#installation)
* [Input File Format](#input-file-format)
* [Usage Examples](#usage-examples)
* [History](#history)
* [Release Notes](#release-notes)
* [License](#license)

## Updates from Previous Releases
We have switched the dependency on the BamTools API to SAMTools, which 
we believe will be more convenient for most users of preseq. Minor bugs 
have been fixed, and algorithms have been refined to more accurately 
construct counts histograms and extrapolate the complexity curve. More 
options have been added to `lc_extrap`. `c_curve` and `lc_extrap` are now both 
under a single binary for easier use, and commands will now be written 
as `preseq lc_extrap [OPTIONS]`. Furthermore, there are updates to the
manual for any minor issues encountered when compiling the preseq
binary.


## Contact Information
Timothy Daley
tdaley@usc.edu
http://smithlabresearch.org

## System Requirements
The preseq software will only run on 64-bit UNIX-like operating 
systems and was developed on Linux systems. The preseq software 
requires a fairly recent C++ compiler (i.e. it must include tr1 
headers). preseq has been compiled and tested on Linux and Mac 
OS X operating systems using GCC v4.1 or greater. 

## Installation
This should be easy: unpack the archive and change into the archive
directory. Then type `make all`. The programs will be in the archive
directory. These can be moved around, and also do not depend on any
dynamic libraries, so they should simply work when executed. If the 
desired input is in .bam format, SAMTools is required. Type 
`make all SAMTOOLS_DIR=/samtools_loc/` to make the programs.

## Input File Format
Input files can be either in BED or BAM file format.  The file should
be sorted by chromosome, start position, strand position, and finally 
strand if in BED format. If the file is in BAM format, then the file
should be sorted using BamTools or SAMTools sort.

## Usage Examples
Each program included in this software package will print a list of
options if executed without any command line arguments. Many of the
programs use similar options (for example, output files are specified
with `-o`). To predict the yield of a future experiment, use `lc_extrap`.
For the most basic usage of `lc_extrap` to compute the expected yield,
use the command:

```
  preseq lc_extrap -o yield_estimates.txt input.bed
```

If the input file is in .bam format, use the command:

```
  preseq lc_extrap -B -o yield_estimates.txt input.bam
```

The yield estimates will appear in yield_estimates.txt, and will be a 
column of future experiment sizes in `TOTAL_READS`, a column of the 
corresponding expected distinct reads in `EXPECTED_DISTINCT`, followed by 
two columns giving the corresponding confidence intervals.  

To investigate the past yield of an experiment, use `c_curve`.  For the
most basic usage, use the command:

```
  preseq c_curve -o estimates.txt input.bed
```

If the input file is in .bam format, use the command:

```
  preseq c_curve -B -o estimates.txt input.bam
```

The estimates will appear in `estimates.txt` with two columns.  The
first column gives the total number of reads in a theoretically
smaller experiment and the second gives the corresponding number of
distinct reads.

## History
preseq was originally developed by Timothy Daley and Andrew Smith 
at the University of Southern California.

## Release Notes

* preseq beta release 1.0.0 (Mar 21, 2014)
  * Included gc_extrap option to predict genomic coverage for
    single cell sequencing experiments.
  * Changed the method of finding optimal continued fraction
    to improve performance for high variable (e.g. RNAseq)
    libraries.
* preseq beta release 0.0.3 (Aug 5, 2013)
  * Prediction of the complexity curve is done using the observed data
    when possible.  Previous versions bootstrapped the histogram and
    used the median estimate.  Bootstrapping is only done to compute
    confidence intervals.
  * Addition of quick mode option with flag `-Q`.  The complexity is 
    predicted with the observed data and bootstrapping is not done,
    speeding up the computation time tremendously.
  * Fixed a bug associated with the `unistd.h` header for GCC versions
    4.7+.  
  * Extensively updated the manual with examples and FAQ


## License
The preseq software for estimating complexity
Copyright &copy; 2013 Timothy Daley and Andrew D Smith and 
the University of Southern California

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
