#!/bin/bash

fnames=(`grep -lR --include '*.log' -e '^TOAST INFO: Workflow completed in' logs/LAT0_CHLAT`)
echo "Moving ${#fnames[@]} successful logs"
for fname in ${fnames[*]}; do
    echo "Moving $fname"
    indir=`dirname $fname`
    infile=`basename $fname`
    outdir=${indir/logs/cleared_logs}
    echo "Moving $fname to $outdir/$infile" 
    mkdir -p $outdir
    mv $fname $outdir
done
echo "Moved ${#fnames[@]} successful logs"

fnames=(`find -L logs -mmin +60 -name '*.log'`)
echo "Removing ${#fnames[@]} stalled logs"
for fname in ${fnames[*]}; do
    echo "Removing $fname"
    indir=`dirname $fname`
    infile=`basename $fname`
    outdir=${indir/logs/failed_logs}
    echo "Moving $fname to $outdir/$infile" 
    mkdir -p $outdir
    mv $fname $outdir
done
echo "Removed ${#fnames[@]} stalled logs"
