#!/bin/bash

#for suffix in "" "_lowcomplexity" "_highcomplexity"; do
for suffix in ""; do
    echo "suffix = ${suffix}"

    fnames=(`grep -lR --include '*.log' -e '^TOAST INFO: Workflow completed in' logs${suffix}/LAT0_CHLAT`)
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

    #fnames=(`grep -lR --include '*.log' -e 'DUE TO TIME LIMIT' logs/LAT0_CHLAT`)
    fnames=(`find -L logs${suffix} -mmin +60 -name '*.log'`)
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
done

exit

fnames=(`grep -LR --include '*.log' -e '^TOAST INFO: Workflow completed in' logs/LAT0_CHLAT`)
echo "Removing ${#fnames[@]} failed logs"
for fname in ${fnames[*]}; do
    echo "Removing $fname"
    indir=`dirname $fname`
    infile=`basename $fname`
    outdir=${indir/logs/failed_logs}
    echo "Moving $fname to $outdir/$infile" 
    mkdir -p $outdir
    mv $fname $outdir
done
echo "Removed ${#fnames[@]} failed logs"
