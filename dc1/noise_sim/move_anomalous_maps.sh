#!/bin/bash

for fname in `awk '{if ($2 > 1000) print $4}' mapstats_f150.txt`; do
    indir=`dirname $fname`
    outdir=${indir/outputs/anomalous\/outputs}
    outdir=`dirname $outdir`
    echo mv $indir $outdir
    #mkdir -p $outdir
    #mv $indir $outdir/
done
