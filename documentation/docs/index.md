# Quick start

You can install preseq as follows using conda:
```console
$ conda install -c bioconda preseq
```
The instructions for installing conda are
[here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). You can obtain preseq source code
[here](https://github.com/smithlabcode/preseq).

This documentation is still under construction... But the PDF
documentation still applies and can be found
[here](https://github.com/smithlabcode/preseq).


# Preseq

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


# Installation

## Using conda

You can install preseq as follows using conda:
```console
$ conda install -c bioconda preseq
```
The instructions for installing conda are
[here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

## From source

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

## Detailed usage

c\_curve
--------

The `c_curve` command is used to compute the expected complexity curve
using a hypergeometric formula (Heck et al., 1975, below). This can be
done with input as mapped reads, a counts file, or a counts histogram
file. The output is a text file with two columns. The first gives the
total number of reads and the second the corresponding number of
distinct reads.

```txt
Heck Jr, K. L., van Belle, G., & Simberloff, D. (1975).
Explicit calculation of the rarefaction diversity measurement
and the determination of sufficient sample size.
Ecology, 56(6), 1459-1461.
```

Command line arguments are as follows:
```txt
  -o, -output   yield output file (default: stdout)
  -s, -step     step size in extrapolations [1000000.000000]
  -v, -verbose  print more information
  -P, -pe       input is paired end read file
  -H, -hist     input is a text file containing the observed histogram
  -V, -vals     input is a text file containing only the observed
                counts
  -B, -bam      input is in BAM format
  -l, -seg_len  maximum segment length when merging paired end bam
                reads [5000]
  -r, -seed     seed for random number generator [408]
```

lc\_extrap
----------

The `lc_extrap` command generates the expected yield for theoretical
larger experiments. It can also gives bounds on the number of distinct
reads in the library and the associated confidence intervals. The
confidence interval is computed by bootstrapping the observed
duplicate counts histogram.  The output is a text file with four
columns:

1. The first is the total number of reads
2. The second is the corresponding average expected number of distinct reads
3. The third is the lower limit of the confidence interval
4. The fourth is the upper limit of the confidence interval

Specifying verbose will print the counts histogram of the input file.

Usage:
```console
$ preseq lc_extrap [OPTIONS] <input-file>
```

Command line arguments are as follows:
```txt
  -o, -output   yield output file (default: stdout)
  -e, -extrap   maximum extrapolation [10000000000.000000]
  -s, -step     extrapolation step size [1000000.000000]
  -n, -boots    number of bootstraps [100]
  -c, -cval     level for confidence intervals [0.950000]
  -x, -terms    maximum terms in estimator [100]
  -v, -verbose  print more info
  -B, -bam      input is in BAM format
  -l, -seg_len  maximum segment length when merging paired end bam
                reads [5000]
  -P, -pe       input is paired end read file
  -V, -vals     input is a text file containing only the observed
                counts
  -H, -hist     input is a text file containing the observed histogram
  -Q, -quick    quick mode (no bootstraps) for confidence intervals
  -D, -defects  no testing for defects
  -r, -seed     seed for random number generator [408]
```

gc\_extrap
----------

The `gc_extrap` command is used to extrapolate the expected number of
bases covered at least once for theoretical larger experiments.  This
is useful for single cell or low input sequencing experiments. Input
format is required to be in mapped read format and we have provided
the tool `to-mr` to convert BAM format files to the mr format. The
output is a text file with four columns.

1. The first is the total number of sequenced and mapped bases
2. The second gives the corresponding expected number of distinct bases covered
3. The third gives the lower limit of the confidence interval
4. The fourth gives the upper limit of the confidence interval

Specifying verbose will print the coverage counts histogram of the
input file.

Usage:
```console
$ preseq gc_extrap [OPTIONS] <input-file>
```

Command line arguments are as follows:
```txt
  -o, -output      coverage yield output file (default: stdout)
  -w, -max_width   max fragment length, set equal to read length for
                   single end reads [10000]
  -b, -bin_size    bin size [10]
  -e, -extrap      maximum extrapolation in base pairs
                   [1000000000000.000000]
  -s, -step        step size in bases between extrapolations
                   [100000000.000000]
  -n, -bootstraps  number of bootstraps [100]
  -c, -cval        level for confidence intervals [0.950000]
  -x, -terms       maximum number of terms [100]
  -v, -verbose     print more information
  -B, -bed         input is in bed format without sequence information
  -Q, -quick       quick mode: run gc_extrap without bootstrapping for
                   confidence intervals
  -D, -defects     defects mode to extrapolate without testing for
                   defects
  -r, -seed        seed for random number generator [408]
```

bound\_pop
----------

The `bound_pop` command is a method for estimating species richness,
the total number of species or classes in the sampled population.
Input format is the same as `lc_extrap`.  Default output is a three
column text file, with the first column containing the estimated
species richness and the second and third containing the estimated
lower and upper confidence intervals. If `bound_pop` is run in quick
mode, then the output is two columns.  The first column will contain
the estimated species richness and the second column will contain the
dimension or order of the approximation.

Usage:
```console
$ preseq bound_pop [OPTIONS] <input-file>
```

Command line arguments are as follows:
```txt
  -o, -output          species richness output file (default: stdout)
  -p, -max_num_points  maximum number of points in quadrature estimates
                       [10]
  -t, -tolerance       numerical tolerance [0.000000]
  -n, -bootstraps      number of bootstraps [500]
  -c, -clevel          level for confidence intervals [0.950000]
  -v, -verbose         print more information
  -P, -pe              input is paired end read file
  -H, -hist            input is a text file containing the observed
                       histogram
  -V, -vals            input is a text file containing only the
                       observed duplicate counts
  -B, -bam             input is in BAM format
  -l, -seg_len         maximum segment length when merging paired end
                       bam reads [5000]
  -Q, -quick           quick mode, estimate without bootstrapping
  -r, -seed            seed for random number generator [408]
```

Examples
========

lc\_extrap
----------

**Using a file of mapped reads**

The mapped reads given as input must be sorted, and can be in BED
format or BAM format if the `-B` flag is specified. We will use the
data file `SRR1041830_5M.bam` provided through a link in the `data`
subdirectory of the preseq source tree. This is a random sample of 5M
reads from a particular sequencing data set.
```console
$ preseq lc_extrap -B -o future_yield.txt SRR1041830_5M.bam
```

Notice that the `-B` flag was specified on the command line. This is
required if the input data is in BAM format. If you run `preseq
lc_extrap` and it does not list the `-B` argument, then your version
of preseq was not built with BAM support.

The output file `future_yield.txt` will contain estimaes of the number
of distinct reads that will be sequenced upon deeper sequencing from
the same library. In this case, we can see some of these predictions
as follows. The first several predictions can be seen as follows:
```console
$ head future_yield.txt
TOTAL_READS EXPECTED_DISTINCT   LOWER_0.95CI    UPPER_0.95CI
0   0   0   0
1000000.0   990651.5    990601.3    990709.5
2000000.0   1963963.0   1963805.8   1964152.2
3000000.0   2921362.5   2921067.3   2921723.9
4000000.0   3863837.0   3863361.6   3864408.8
5000000.0   4792170.5   4791451.2   4792972.3
6000000.0   5706990.8   5705989.7   5708096.4
7000000.0   6608846.9   6607533.2   6610298.1
8000000.0   7498214.9   7496569.6   7500056.3
```

Note that the first line above gives column headings. The second line
corresponds to sequencing zero reads (and is as expected). The
remaining lines give estimates for what we would expect if we
sequenced the corresponding number of reads from the library. If we
sequenced 8M reads, then we expect 7498214.9 of them to be
unique. Here the confidence intervals are quite narrow, but that's not
always the case, especially in very low-complexity libraries.

We can examine the estimates for sequencing 50M, 500M and 5G reads as
follows:
```console
$ grep ^50000000 future_yield.txt
50000000.0      36786623.6  35117618.0  37112536.8
500000000.0     116108773.9 97165891.1  137121043.6
5000000000.0    148401854.8 118080946.1 198715938.8
```

Above the second column shows the expected distinct reads for
sequencing 10x, 100x and 1000x deeper from the same library. In this
case, we still sequence new unique reads after 500M, but but the
fraction unique at 5G is under 3%. Looking more closely around the 5G
estimate:
```console
$ grep -2 ^5000000000 future_yield.txt
4998000000.0    148400017.5 118079815.0 198711828.1
4999000000.0    148400936.4 118080380.7 198713883.8
5000000000.0    148401854.8 118080946.1 198715938.8
5001000000.0    148402772.9 118081511.4 198717993.0
5002000000.0    148403690.7 118082076.4 198720046.3
```
we see that the library likely only contains 150M unique molecules.

This example uses a sorted read mapped reads file as the data for the
initial experiment.  The default step size between yield estimates is
1M, the default confidence interval level is 95%, and the default
total extrapolation is to 10 billion (10G) sequenced reads.

**Using a mapped reads file and the verbose option**

```console
$ preseq lc_extrap -v -B -o future_yield.txt SRR1041830_5M.bam
```

As preseq is running, information will print to the terminal screen
(standard error) that includes a read counts histogram of the input
file. This truncates after the first value that has zero
observations. Included here is the first 20 lines of what would be
observed:
```txt
BAM_INPUT
TOTAL READS     = 5083564
DISTINCT READS  = 4.86904e+06
DISTINCT COUNTS = 19
MAX COUNT       = 19
COUNTS OF 1     = 4.67849e+06
MAX TERMS       = 18
OBSERVED COUNTS (20)
1   4678490
2   173897
3   12770
4   2345
5   772
6   355
7   156
8   82
9   55
10  37
11  24
12  21
```

** Using a sorted mapped reads with options**

```console
$ preseq lc_extrap -B -e 15000000 -s 500000 -b 90 -c .90 -o future_yield.txt SRR1041830_5M.bam
```

This will use a step size of 500k when making estimates, and the
estimates ending at 15M. The confidence intervals are now at a level
of 90%. This can be seen as follows:
```console
$ head future_yield.txt
TOTAL_READS EXPECTED_DISTINCT   LOWER_0.9CI UPPER_0.9CI
0   0   0   0
500000.0    497614.0    497602.0    497630.0
1000000.0   990651.5    990609.0    990703.1
1500000.0   1479369.5   1479288.0   1479477.3
2000000.0   1963963.0   1963834.6   1964139.1
2500000.0   2444581.0   2444401.0   2444831.2
3000000.0   2921362.5   2921116.3   2921682.2
3500000.0   3394409.5   3394093.8   3394798.0
4000000.0   3863837.0   3863430.4   3864324.8
$ tail -5 future_yield.txt
12500000.0  11356784.0  11353132.2  11360765.0
13000000.0  11771818.5  11767804.8  11776196.8
13500000.0  12184179.7  12179952.7  12188987.7
14000000.0  12594028.0  12589384.9  12599334.0
14500000.0  13001405.8  12996034.3  13007337.6
```

**Using a histogram or read counts as input**

`lc_extrap` allows the input file to be an observed histogram. An
example of the format of this histogram is (truncated at 20 lines):
```txt
1       47879524
2       16751354
3       7193291
4       3552776
5       1935072
6       1120980
7       682786
8       430623
9       280089
10      186473
11      127529
12      89218
13      63197
14      45779
15      33644
16      25155
17      18760
18      14371
19      11048
20      8667
```

The following command will give the same output as the first example
above (using default parameters):
```console
$ preseq lc_extrap -o future_yield.txt -H histogram.txt
```
Assuming that `histogram.txt` has the proper format and includes the
frequencies of read counts from `SRR1041830_5M.bam`.

Similarly, both `lc_extrap` and `c_curve` allow the option to input
read counts (text file should contain ONLY the observed counts in a
single column). For example, if a dataset had the following counts
histogram:
```txt
1      4
2      3
3      1
```
then the corresponding input file of just read counts might be
```txt
1
2
1
1
3
2
2
1
```

This particular example involves counts that are far to small for
preseq to make any estimates, but they should illustrate the
relationship between a file of counts, and the counts histogram
format.

When using these count "values" the command should be run with the `-V`
flag (not to be confused with the lowercase `-v` for verbose mode):
```console
$ preseq lc_extrap -o future_yield.txt -V counts.txt
```

gc\_extrap
----------

The `gc_extrap` command is designed for estimating genome coverage
from more deeply sequencing a library, particularly in single cell
whole genome sequencing experiments. For illustrative purposes we will
examine an MDA (multiple displacement amplification) whole genome
sequencing experiment, SRA accession SRR1777281. This experiment has
5.76 million paired end 101 base reads. We mapped the experiment with
bowtie2 v0.0-beta7 under default parameters. This resulted in 3.63
million concordantly mapped fragment pairs and 2.1 million
disconcordantly mapped fragments.

The first step is to convert the sorted bam file to mr format
and sort it.
```console
$ to-mr -o SRR1777281_bwt2.mr -L 10000 SRR1777281_bwt2.sort.bam
$ sort -k 1,1 -k 2,2n -k 3,3n SRR1777281_bwt2.mr > SRR1777281_bwt2.sort.mr
```
The resulting mapped reads file has 813 million bases total (note that
bases covered by two fragments of the same read are only counted once)
and 410 million covered bases in the genome.

By default, `gc_extrap` divides the genome into 10 base pair non-overlapping bins.
In default mode, the running time of `gc_extrap` was under 12 minutes.
```txt
LOADING READS
MAPPED READ FORMAT
TOTAL READS         = 5726883
BASE STEP SIZE      = 1e+08
BIN STEP SIZE       = 1e+07
TOTAL BINS          = 8.1325e+07
BINS PER READ       = 14.2006
DISTINCT BINS       = 4.09582e+07
TOTAL BASES         = 8.1325e+08
TOTAL COVERED BASES = 4.09582e+08
MAX COVERAGE COUNT  = 79775
COUNTS OF 1         = 3.13723e+07
OBSERVED BIN COUNTS (79776)
1       3.13723e+07
2       6.97799e+06
3       1.72101e+06
```

The first several lines of the output file were as follows:
```txt
TOTAL_BASES    EXPECTED_COVERED_BASES    LOWER_95%CI    UPPER_95%CI
0    0     0     0
100000000.0     64522380.0     63075879.1     66002053.1
200000000.0     123422455.0     120747705.7     126156454.1
300000000.0     178054120.0     174319619.1     181868626.2
400000000.0     229008295.0     224352188.5     233761032.3
500000000.0     276727080.0     271265011.6     282299130.1
```
The final line of the output was:
```txt
999900000000.0     1826891418.6     1621208958.1     2058668772.4
```

To run `gc_extrap` at single base resolution, the option
`-b 1` is required.  This results in a significant increase
in the running time of the algorithm.  For this case the
running time was 113 minutes. The following is the command
along with the first several lines of information reported
to the terminal:
```console
$ preseq gc_extrap SRR1777281_bwt2.sort.mr -o SRR1777281_bwt2_1bp_gc_extrap.txt -b 1 -v
LOADING READS
MAPPED READ FORMAT
TOTAL READS         = 5726883
BASE STEP SIZE      = 1e+08
BIN STEP SIZE       = 1e+08
TOTAL BINS          = 8.13236e+08
BINS PER READ       = 142.003
DISTINCT BINS       = 4.09454e+08
TOTAL BASES         = 8.13236e+08
TOTAL COVERED BASES = 4.09454e+08
MAX COVERAGE COUNT  = 80028
COUNTS OF 1         = 3.13527e+08
OBSERVED BIN COUNTS (80029)
1       3.13527e+08
2       6.98288e+07
3       1.722e+07
```

The output looks like this:
```txt
TOTAL_BASES    EXPECTED_COVERED_BASES    LOWER_95%CI    UPPER_95%CI
0    0     0     0
100000000.0    64680427.0    64284593.4    65078698.0
200000000.0    123709880.0    122978429.6    124445680.9
300000000.0    178447901.0    177427065.8    179474609.6
400000000.0    229488768.0    228216803.2    230767822.1
500000000.0    277279369.5    275788078.2    278778724.8
```
with the final line:
```txt
999900000000.0    1838021604.0    1682515315.7    2007900543.3
```
