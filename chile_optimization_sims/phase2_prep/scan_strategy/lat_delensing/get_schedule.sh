#/bin/bash

toast_ground_schedule \
    @schedule_lat_delensing.par \
    @patches_lat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --sun-avoidance-angle 30 \
    --moon-avoidance-angle 0 \
    --out schedule_lat_delensing.split1.txt \
    >& schedule_lat_delensing.split1.log &

toast_ground_schedule \
    @schedule_lat_delensing.par \
    @../sat/patches_sat.txt  \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --sun-avoidance-angle 30 \
    --moon-avoidance-angle 0 \
    --out schedule_lat_delensing.split2.txt \
    >& schedule_lat_delensing.split2.log &

#     --patch North_Direction,SIDEREAL,100.0,330.00,30.00,67.00,130,160,15
