#/bin/bash

#export TOAST_LOGLEVEL=DEBUG


python make_targets.py

toast_ground_schedule \
    @schedule_lat_delensing.par \
    --start "2030-01-01 00:00:00" \
    --stop "2031-01-01 00:00:00" \
    @patches.txt \
    --time-step 1200 \
    --out schedule_lat_delensing_max.txt \
    >& schedule_lat_delensing_max.log &

# python gapfill_schedule.py schedule_sat.sun90max.txt schedule_sat.sun45max.txt schedule_sat.sun45_supplement.txt
