# preseq

Under construction... But the PDF documentation still applies and can
be found [here](https://github.com/smithlabcode/preseq).

The preseq package is aimed to help researchers design and optimize
sequencing experiments by using population sampling models to infer
properties of the population or the behavior under deeper sampling
based upon a small initial sequencing experiment.  The estimates can
then be used to examine the utility of further sequencing, optimize
the sequencing depth, or to screen multiple libraries to avoid low
complexity samples.

The four main commands are `c_curve`, `lc_extrap`, `gc_extrap`, and
`bound_pop`. The `c_curve` command interpolates the expected
complexity curve based upon a hypergeometric formula and is primarily
used to check predictions from `lc_extrap` and `gc_extrap`. It will
show you what your data looks like (as a curve) within the range of
data you already have. The `lc_extrap` command uses rational function
approximations to Good & Toulmin's (1956) non-parametric empirical
Bayes estimator to predict the library complexity of future
experiments, in essence looking into the future for hypothetical
experiments. This is the paper:

```txt
Good IJ & Toulmin GH (1956).
The number of new species, and the increase in population coverage,
when a sample is increased.
Biometrika, 43(1-2),45-63.
```

The `gc_extrap` command uses a similar approach as `lc_extrap` to
predict the genome coverage, i.e. the number of bases covered at least
once, from deeper sequencing in a single cell or low input sequencing
experiment based on the observed coverage counts.  An option is
available to predict the coverage based on binned coverage counts to
speed up the estimates. `gc_extrap` requires mapped read or bed format
input, so the tool `to-mr` is provided to convert bam format read to
mapped read format.

`bound_pop` uses a non-parametric moment-based approach to
conservatively estimate the total number of classes in the sample,
also called the species richness of the population that is sampled.


## Installation

**Download**

You can obtain preseq [here](https://github.com/smithlabcode/preseq).
Only clone the repo if you are planning to modify or reuse the code,
or have some other interesting plans and know what you are doing. Most
users should download a release. The `README.md` file in the GitHub
repo explains how to download a release, and not simply the source
code files.

Preseq runs on Unix-based systems, notably macOS and Linux. If the
input file is in BAM format, the
[HTSLib](http://www.htslib.org/download/) library is required. If the
input is a text file of counts in a single column, a "counts
histogram" or is in BED format, then HTSLib is not required. Preseq
has been tested on Linux and Mac OS-X.

***Installation**

Instructions on how to install preseq are included in the `README.md`
file in root of the source tree for preseq.  If you have problems
installing preseq, please post an issue on the GitHub page or contact
Andrew Smith <andrewds@usc.edu>.

Using preseq
============

Basic usage
-----------

To generate the complexity curve of a genomic library from a read file
in BED or BAM format or a duplicate count file, use the `c_curve`
command. Use `-o` to specify the output file name.
```console
$ preseq c_curve -o complexity_output.txt input.bed
```

In the command above (and everywhere in these docs) we assume that the
`preseq` program is in your path

To predict the complexity curve of a sequencing library using an
initial experiment in BED format, use the `lc_extrap` command. The
required options are `-o` to specify the output of the yield estimates
and the input file, which is either a BED file sorted by chromosome,
start position, end position, and strand or a BAM file sorted with the
samtools sort function. Additional options are available and are
detailed below.
```console
$ preseq lc_extrap -o future_yield.txt input.bed
```

For a low input sequencing experiment the genomic coverage is highly
variable and uncertain function of sequencing depth.  Some regions may
be missing due to locus dropout or preferentially amplified during
whole genome amplification. The `gc_extrap` command allows the level
of genomic coverage from deep sequencing to be predicted based on an
initial sample.  The input file format need to be a mapped read
(`.mr`) or BED, sorted by chromosome, start position, and end then
position.  Additional options are available and are detailed below.
```console
$ preseq gc_extrap -o future_coverage.txt input.mr
```

File Format
===========

**Sorted read files in BED or BAM format**

Input files are sorted mapped read files in BED or BAM format, or a
text file consisting of one column giving the observed read counts.
The programs require that BED files are sorted by chromosome, start
position, end position, and strand.  This can be achieved by using the
command line function sort as follows:
```console
$ sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 input.bed > input.sort.bed
```

BAM format read files should be sorted by chromosome and start
position. This can be done with the SAMTools sort function. If the
input is in BAM format, then the flag `-B` must be included.

If the input is for paired end mapped reads, the option `-P` can be
set.  In this case concordantly mapped reads and disconcordantly
mapped fragments are counted.  This means that both ends of a
disconcordantly mapped read will each be counted separately.  If a
large number of reads are disconcordant, then the default single end
should be used or the disconcordantly mapped reads removed prior to
running `preseq`.  In this case only the mapping location of the first
mate is used as the unique molecular
identifier. See this paper for details:
```txt
Kivioja T, Vähärautio A, Karlsson K, Bonke M, Enge M, Linnarsson S & Taipale J (2012).
Counting absolute numbers of molecules using unique molecular identifiers.
Nature Methods, 9(1),72-74.
```

**Text files of observed read counts**

For more general applications `preseq` allows the input
to be a text file of observed read counts, one count per
line. To specify this input, the option `-V` must be set.

Such a text file can typically be constructed using shell commands or
very simple scripts. Take for example an unmapped sequencing
experiment in FASTQ format. To predict the complexity, the unique
molecular identifier needs to use only the observed sequence.  For
instance, a unique molecular identifier used may be the first 20 bases
in the observed sequence. A command line script to construct the
counts would then be
```console
$ awk '{if (NR%4==2) print substr($0,1,20);}' input.fastq | sort | \
    uniq -c | awk '{print $1}' > counts.txt
```
More complicated unique molecular identifiers can be used, such as
mapping position plus a random barcode, but are too complicated to
explain here. For questions with such usage, please contact
us at <andrewds@usc.edu>

**Mapped read format for `gc_extrap`**

`gc_extrap` does not allow for input files to be in BAM format.  We
have found that some mappers give inconsistent SAM flags for paired
end reads, preventing efficient merging of reads in the proper order.
We provide the tool `to-mr` to convert SAM or BAM format files to
MR format.  The MR or BED format file needs to be sorted by
chromosome, start, and end position before input into `gc_extrap`.
