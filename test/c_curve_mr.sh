#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0

prog=./preseq
infile=data/SRR1003759_5M_subset.mr
outdir=c_curve_mr_out
if [[ -e "${infile}" ]]; then
    mkdir -p ${outdir}
    ${prog} c_curve --step 1000 -o ${outdir}/out.txt ${infile}
    x=$(md5sum --ignore-missing -c test/md5sum.txt | \
            grep "${outdir}" | \
            grep -c "OK$")
    if [[ "${x}" != "1" ]]; then
        exit 1;
    fi
    rm -r ${outdir}
else
    echo "${infile} not found; skipping remaining tests";
    exit 77;
fi
