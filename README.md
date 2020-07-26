This is the README file for the preseq package. The preseq package is
aimed at predicting the yield of distinct reads from a genomic library
from an initial sequencing experiment. The estimates can then be used
to examine the utility of further sequencing, optimize the sequencing
depth, or to screen multiple libraries to avoid low complexity
samples.

SYSTEM REQUIREMENTS:
========================================================================
The preseq software will only run on 64-bit UNIX-like operating
systems and was developed on both Linux and Mac. The preseq software
requires a C++ compiler that supports C++11. The GNU Scientific Library
is also required. It can be installed using `apt` on Linux, using `brew` 
on macOS, or from source available [here](http://www.gnu.org/software/gsl).

INSTALLATION:
========================================================================
### Installing from a release

1. Download preseq-3.0.tar.gz from the releases tab of this repository.
2. Unpack the archive:
```
$ tar -zxvf preseq-3.0.tar.gz
```
3. Move into the preseq directory and create a build directory:
```
$ cd preseq-3.0
$ mkdir build && cd build
```
4. Run the configuration script:
```
$ ../configure
```
If you do not want to install preseq system-wide, or if you do
not have admin privileges, specify a prefix directory:
```
$ ../configure --prefix=/some/reasonable/place
```
Finally, if you want to build with HTSlib support (for the `to-mr`
program) then you need to specify the following:
```
$ ../configure --enable-hts
```
And if you installed HTSlib yourself in some non-standard directory,
you must specify the location like this:
```
$ ../configure --enable-hts CPPFLAGS='-I /path/to/htslib/headers' \
    LDFLAGS='-L/path/to/htslib/lib'
```
5. Compile and install the tools:
```
$ make
$ make install
```

### Installing from source

Developers looking to use the latest commits can compile the cloned
repository using the `Makefile` within the `src` directory. The
process is simple:
```
$ cd src/
$ make
```
If the desired input is in `.bam` format, `htslib` is required. Type
```
make HAVE_HTSLIB=1 all
```
to make the programs. The HTSLib library can be obtained here:
http://www.htslib.org/download

INPUT FILE FORMATS:
========================================================================
The input to preseq can be in 3 general formats:
1. Mapped read locations in BED or BAM file format. The file should be
   sorted by chromosome, start position, end position, and finally
   strand if in BED format. If the file is in BAM format, then the
   file should be sorted using `bamtools` or `samtools sort`.
2. The "counts histogram" which will have, for each count 1,2,..., the
   number of unique "species" (e.g. reads, or anything else) that
   appear with that count. Examples can be found in the data directory
   within the preseqR subdirectory. Note these should not have a count
   for "0", and they should not have any header above the counts. Just
   two columns of numbers, with the first column sorted and unique.
3. The counts themselves, so just a file with one count on each
   line. These will be made into the "counts histogram" inside preseq
   right away.

USAGE EXAMPLES:
========================================================================
Each program included in this software package will print a list of
options if executed without any command line arguments. Many of the
programs use similar options (for example, output files are specified
with '-o'). To predict the yield of a future experiment, use `lc_extrap`.
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
corresponding expected distinct reads in `EXPECTED_DISTINCT`, followed
by two columns giving the corresponding confidence intervals.

To investigate the past yield of an experiment, use `c_curve`.  For the
most basic usage, use the command:
```
preseq c_curve -o estimates.txt input.bed
```
If the input file is in .bam format, use the command:
```
preseq c_curve -B -o estimates.txt input.bam
```
The estimates will appear in estimates.txt with two columns.  The
first column gives the total number of reads in a theoretically
smaller experiment and the second gives the corresponding number of
distinct reads.

UPDATES TO VERSION 3.0
========================================================================
The main change to this version is that if BAM/SAM format will be used
as input, the HTSLib library must be installed on the system when
preseq is built. Installation instructions below have been updated
correspondingly. We also updated to use C++11, so a more recent
compiler is required, but these days C++11 is usually supported.

UPDATES TO VERSION 2.0.3
========================================================================
A bug in defect mode was fixed and a rng seed was added to allow for
reproducibility.

UPDATES TO VERSION 2.0.0
========================================================================
We have added a new module, `bound_pop`, to estimate a lower bound of
the population sampled from.  Interpolation is calculated by
expectation rather than subsampling, dramatically improving the speed.

UPDATES TO VERSION 1.0.2
========================================================================
We have switched the dependency on the BamTools API to SAMTools, which
we believe will be more convenient for most users of preseq. Minor
bugs have been fixed, and algorithms have been refined to more
accurately construct counts histograms and extrapolate the complexity
curve. More options have been added to `lc_extrap`. `c_curve` and
`lc_extrap` are now both under a single binary for easier use, and
commands will now be written as `preseq lc_extrap [OPTIONS]`
Furthermore, there are updates to the manual for any minor issues
encountered when compiling the preseq binary.

We released an R package called
[preseqR](http://cran.r-project.org/web/packages/preseqR/index.html)
along with preseq. This makes most of the preseq functionality
available in the R statistical environment, and includes some new
functionality. The preseqR directory contains all required source code
to build this R package.

CONTACT INFORMATION:
========================================================================
Andrew D. Smith
andrewds@usc.edu

Timothy Daley
tdaley@stanford.edu

http://smithlabresearch.org

HISTORY
========================================================================
Preseq was originally developed by Timothy Daley and Andrew D. Smith
at University of Southern California.

LICENSE
========================================================================
The preseq software for estimating complexity Copyright (C) 2014-2020
Timothy Daley and Andrew D Smith and Chao Deng and the University of
Southern California

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
