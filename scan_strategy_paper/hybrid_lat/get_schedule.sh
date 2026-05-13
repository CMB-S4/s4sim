#/bin/bash

#export TOAST_LOGLEVEL=DEBUG

if [[ ! -e schedule_lat_deep.txt ]]; then

    python make_targets.py

    toast_ground_schedule \
        @schedule_lat_delensing.par \
        @patches.txt \
        --time-step 1200 \
        --out schedule_lat_deep.txt \
        >& schedule_lat_deep.log &
fi

if [[ ! -e schedule_lat_wide.txt ]]; then
    toast_ground_schedule \
        @schedule_lat_wide.par \
        --out schedule_lat_wide.txt \
        >& schedule_lat_wide.log &
fi

wait

toast_gapfill_schedule schedule_lat_deep.txt schedule_lat_wide.txt schedule_lat_supplement.txt --supplement_only
