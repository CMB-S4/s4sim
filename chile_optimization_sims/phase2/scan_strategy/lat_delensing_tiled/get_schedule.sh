#/bin/bash

toast_ground_schedule \
    @schedule_lat_delensing.par \
    @patches.txt \
    --out schedule_lat_delensing_tiled.txt \
    >& schedule_lat_delensing_tiled.log &

#    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
