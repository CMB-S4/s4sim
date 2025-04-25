#/bin/bash

# The new schedule observes at el=35 rather than el=40.
# The throw is 140 degrees rather than 120 degrees which
# increases the scan rate modulation factor from 2 to 2.92.
# An extra 10 degrees in throw would make the factor 3.86.

toast_ground_schedule \
    @schedule_lat_wide.par \
    --out schedule_lat_wide.txt \
    >& schedule_lat_wide.log &

# --block-out 01/01-04/01
# python gapfill_schedule.py schedule_lat_wide.txt ../lat_delensing_max/schedule_lat_delensing_max.txt schedule_lat_delensing_supplement.txt
