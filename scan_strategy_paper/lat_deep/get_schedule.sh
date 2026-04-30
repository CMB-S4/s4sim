#/bin/bash

#export TOAST_LOGLEVEL=DEBUG


python make_targets.py

toast_ground_schedule \
    @schedule_lat_delensing.par \
    @patches.txt \
    --time-step 1200 \
    --out schedule_lat_deep.txt \
    >& schedule_lat_deep.log &
