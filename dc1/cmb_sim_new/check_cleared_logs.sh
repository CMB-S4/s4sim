#!/bin/bash

#for suffix in "" "_lowcomplexity" "_highcomplexity"; do
for suffix in "_highcomplexity"; do
    echo "suffix = ${suffix}"

    fnames=(`grep -LR --include '*.log' -e '^TOAST INFO: Workflow completed in' cleared_logs${suffix}/LAT0_CHLAT`)
    echo "Removing ${#fnames[@]} failed logs"
    for fname in ${fnames[*]}; do
        echo "Removing $fname"
        indir=`dirname $fname`
        infile=`basename $fname`
        outdir=${indir/cleared_logs/failed_logs}
        echo "Moving $fname to $outdir/$infile" 
        mkdir -p $outdir
        mv $fname $outdir
    done
    echo "Removed ${#fnames[@]} failed logs"
done
