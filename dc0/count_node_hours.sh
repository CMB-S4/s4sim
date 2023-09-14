#!/bin/bash

workdir=$PWD

cd cleared_logs
for TELE in *; do
    echo "$(date) : $TELE"
    cd $TELE
    for band in *; do
        echo "$(date) : $band"
        cd $band
        fname_out=${workdir}/times_${TELE}_${band}.txt
        rm -f $fname_out
        echo "$(date) : Writing $fname_out"
        touch $fname_out
        for fname in ${workdir}/cleared_logs/${TELE}/${band}/*.log; do
            # Find the first instance of the search string
            nnode=`grep -m 1 "^TOAST INFO: Executing workflow" $fname | awk '{print $6 * $11 / 256}'`
            if [[ "$nnode" = "" ]]; then
                echo "Failed to count nodes in $fname"
                continue
            fi
            # Search for the complete time in the last 10 lines
            elapsed=`tail $fname | grep -m 1 "^TOAST INFO: Workflow completed in" | awk '{print $6}'`
            if [[ "$elapsed" = "" ]]; then
                # Search for the complete time in the entire log
                elapsed=`grep -m 1 "^TOAST INFO: Workflow completed in" $fname | awk '{print $6}'`
                if [[ "$elapsed" = "" ]]; then
                    # Can't do it. Give up on this log
                    echo "Failed to get elapsed time in $fname"
                    continue
                fi
            fi
            echo "$fname $nnode $elapsed" >> $fname_out
        done
        cd ..
    done
    cd ..
done

exit

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
    for fname in staged_logs/LAT0_CHLAT/$band/*.log; do
        nnode=`head $fname | grep "^TOAST INFO: Executing workflow" | awk '{print $6 * $11 / 256}'`
        elapsed=`tail $fname | grep "^TOAST INFO: Workflow completed" | awk '{print $6}'`
        echo "$fname $nnode $elapsed" >> $fname_out
    done
done
