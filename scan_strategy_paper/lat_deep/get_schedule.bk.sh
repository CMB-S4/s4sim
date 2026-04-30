#/bin/bash

#export TOAST_LOGLEVEL=DEBUG


python make_targets.bk.py

toast_ground_schedule \
    @schedule_lat_delensing.par \
    @patches.bk.txt \
    --time-step 1200 \
    --out schedule_lat_deep.bk.txt \
    >& schedule_lat_deep.bk.log &
