#!/bin/bash

fnames=(`grep -lR --include '*.log' -e '^TOAST INFO: Workflow completed in' logs/LAT0_CHLAT`)
echo "Moving ${#fnames[@]} successful logs"
for fname in ${fnames[*]}; do
    echo "Moving $fname"
    indir=`dirname $fname`
    outdir=${indir/logs/cleared_logs}
    echo "Moving $fname to $outdir" 
    mkdir -p $outdir
    mv $fname $outdir
done

fnames=(`grep -LR --include '*.log' -e '^TOAST INFO: Workflow completed in' logs/LAT0_CHLAT`)
echo "Removing ${#fnames[@]} failed logs"
for fname in ${fnames[*]}; do
    echo "Removing $fname"
    indir=`dirname $fname`
    outdir=${indir/logs/failed_logs}
    echo "Moving $fname to $outdir" 
    mkdir -p $outdir
    mv $fname $outdir
done
