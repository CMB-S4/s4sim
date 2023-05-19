#!/bin/bash

for band in f030 f040 f090 f150 f220 f280; do
    echo $band
    fname_out=times_${band}.txt
    rm -f $fname_out
    touch $fname_out
    for fname in cleared_logs/LAT0_CHLAT/$band/*.log; do
        nnode=`head $fname | grep "^TOAST INFO: Executing workflow" | awk '{print $6 * $11 / 256}'`
        elapsed=`tail $fname | grep "^TOAST INFO: Workflow completed" | awk '{print $6}'`
        echo "$fname $nnode $elapsed" >> $fname_out
    done
done
