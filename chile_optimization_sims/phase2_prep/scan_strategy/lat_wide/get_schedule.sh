#/bin/bash

toast_ground_schedule \
    @schedule_lat_wide.par \
    --sun-avoidance-angle 30 \
    --moon-avoidance-angle 0 \
    --out schedule_lat_wide.split1.txt \
    >& schedule_lat_wide.split1.log
