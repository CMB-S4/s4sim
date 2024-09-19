#/bin/bash

toast_ground_schedule \
    @schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --patch North_Direction,SIDEREAL,100.0,330.00,30.00,67.00,130,160,15 \
    --out schedule_sat.txt \
    >& schedule_sat.log
