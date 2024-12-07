#/bin/bash

toast_ground_schedule \
    @schedule_lat_wide.par \
    --out schedule_lat_wide.txt \
    >& schedule_lat_wide.log &

# --block-out 01/01-04/01
