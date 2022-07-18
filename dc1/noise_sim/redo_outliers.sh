#!/bin/bash

#for band in f220 f280; do
for band in f090 f150; do
    for bad in `awk '{if ($2 > 20000) print $1}' maps_${band}.txt.stats`; do
        echo "Moving $bad $band to failed logs"
        mv cleared_logs/LAT0_CHLAT/${band}/${bad}.log failed_logs/LAT0_CHLAT/${band}/
    done
done
