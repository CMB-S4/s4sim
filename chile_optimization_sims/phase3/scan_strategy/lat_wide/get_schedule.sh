#/bin/bash

# The new schedule observes at el=35 rather than el=40.
# The throw is 140 degrees rather than 120 degrees which
# increases the scan rate modulation factor from 2 to 2.92.
# An extra 10 degrees in throw would make the factor 3.86.

toast_ground_schedule \
    @schedule_lat_wide.par \
    --out schedule_lat_wide.txt \
    >& schedule_lat_wide.log &

python gapfill_schedule.py ../lat_delensing_sun90bk/schedule_lat_delensing_sun90bk.txt schedule_lat_wide.txt schedule_lat_wide_supplement.txt
