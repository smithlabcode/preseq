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
requires a C++ compiler that supports C++11.

INSTALLATION:
========================================================================
### Installing from a release

1. Download `preseq-x.tar.gz` from the releases tab of this repository.
2. Unpack the archive:
```
$ tar -zxvf preseq-x.tar.gz
```
3. Move into the preseq directory and create a build directory:
```
$ cd preseq-x
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
The HTSLib library can be obtained here:
http://www.htslib.org/download.

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
with '-o').

We have provided a data directory to test each of our programs.
Change to the `data` directory and try some of our commands.
To predict the yield of a future experiment, use `lc_extrap`.
For the most basic usage of `lc_extrap` to compute the expected yield,
use the command on the following data:
```
preseq lc_extrap -o yield_estimates.txt SRR1003759_5M_subset.mr
```
If the input file is in `.bam` format, use the `-B` flag:
```
preseq lc_extrap -B -o yield_estimates.txt SRR1106616_5M_subset.bam
```
For the counts histogram format, use the `-H` flag:
```
preseq lc_extrap -H -o yield_estimates.txt SRR1301329_1M_read.txt
```

The yield estimates will appear in yield_estimates.txt, and will be a
column of future experiment sizes in `TOTAL_READS`, a column of the
corresponding expected distinct reads in `EXPECTED_DISTINCT`, followed
by two columns giving the corresponding confidence intervals.

To investigate the past yield of an experiment, use `c_curve`.
`c_curve` can take in the same file formats as `lc_extrap` by using
the same flags. The estimates will appear in estimates.txt with two
columns.  The first column gives the total number of reads in a
theoretically smaller experiment and the second gives the
corresponding number of distinct reads.

`bound_pop` provides an estimate for the species richness of the
sampled population. The input file formats and corresponding flags are
identical to `c_curve` and `lc_extrap`. The output provides the median
species richness in the first column and the confidence intervals in
the next two columns.

Finally, `gc_extrap` predicts the expected genomic coverage for a
future experiment.  It produces the coverage in an output format
identical to `lc_extrap`. `gc_extrap` can only take in files in BED
and mapped reads format (using the `-B` flag for BED):
```
preseq gc_extrap -B -o coverage_estimates.txt SRR1003759_5M_subset.mr
```

More data is available in the `additional_data.txt` file in the `data`
directory.  For an extended write-up on our programs, please read the
manual in the `docs` directory.

UPDATES TO VERSION 3.1.0
========================================================================
A mode `pop_size` has been added that uses the continued fraction
approximation to the Good-Toulmin model and extrapolates as far as
possible. Although `bound_pop` provides a good and reliable
lower-bound, this new mode will give a more accurate estimate of the
population size (e.g. total number of distinct molecules). It's not
perfect yet, and in some cases if the population is more than a
billion times larger than the sample, it will still only give a lower
bound. But it works well on most data sets.

UPDATES TO VERSION 3.0.2
========================================================================
GSL has been completely removed, and a data directory has been added
for users to test our programs.

UPDATES TO VERSION 3.0.1
========================================================================
We no longer require users to have GSL for all modules except for
`bound_pop`. Users interested in using `bound_pop` can install GSL and
follow the instructions above to configure with GSL.

UPDATES TO VERSION 3.0
========================================================================
The main change to this version is that if BAM/SAM format will be used
as input, the HTSLib library must be installed on the system when
preseq is built. Installation instructions above have been updated
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
