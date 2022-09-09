#/bin/bash

#python3 make_lat_tiles.py

export TOAST_LOGLEVEL=DEBUG

# Run all the schedules

# --block-out 01/15-03/15

toast_ground_schedule \
    @schedule_lat.par \
    @patches_lat.txt \
    >& get_schedule.log &
