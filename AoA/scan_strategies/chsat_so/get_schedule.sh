#/bin/bash

python3 make_sat_tiles.py

export TOAST_LOGLEVEL=DEBUG

# Run all the schedules

# --block-out 01/15-03/15

toast_ground_schedule \
    @schedule_sat.par \
    @patches_sat.txt \
    >& get_schedule.log &
