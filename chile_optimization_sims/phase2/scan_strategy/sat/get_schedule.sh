#/bin/bash

#toast_ground_schedule \
#    @schedule_sat.par \
#    --block-out "01/01-04/01" \
#    @patches_sat.txt \
#    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
#    --out schedule_sat.w_break.txt \
#    >& schedule_sat.w_break.log &

toast_ground_schedule \
    @schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --sun-avoidance-angle 90 \
    --out schedule_sat.sun90.txt \
    >& schedule_sat.sun90.log &

toast_ground_schedule \
    @schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --sun-avoidance-angle 45 \
    --out schedule_sat.sun45.txt \
    >& schedule_sat.sun45.log &
