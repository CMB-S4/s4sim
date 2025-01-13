#/bin/bash

toast_ground_schedule \
    @pole_schedule_sat.par \
    --sun-avoidance-angle 90 \
    --out schedule_sat.sun90.txt \
    >& schedule_sat.sun90.log &

toast_ground_schedule \
    @pole_schedule_sat.par \
    --sun-avoidance-angle 45 \
    --out schedule_sat.sun45.txt \
    >& schedule_sat.sun45.log &
