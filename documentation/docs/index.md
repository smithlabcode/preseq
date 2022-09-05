# preseq

Under construction... But the PDF documentation still applies and can
be found [here](https://github.com/smithlabcode/preseq).

The `preseq` package is aimed to help researchers design and optimize
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
approximations to Good \& Toulmin's (1956) non-parametric empirical
Bayes estimator to predict the library complexity of future
experiments, in essence looking into the future for hypothetical
experiments.

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
