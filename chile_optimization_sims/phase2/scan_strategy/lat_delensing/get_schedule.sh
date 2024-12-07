#/bin/bash

toast_ground_schedule \
    @schedule_lat_delensing.par \
    @patches_lat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --out schedule_lat_delensing.txt \
    >& schedule_lat_delensing.log &

#toast_ground_schedule \
#    @schedule_lat_delensing.par \
#    @../sat/patches_sat.txt  \
#    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
#    --out schedule_lat_delensing_core.txt \
#    >& schedule_lat_delensing_core.log

#     --patch North_Direction,SIDEREAL,100.0,330.00,30.00,67.00,130,160,15

# --block-out 01/01-04/01
