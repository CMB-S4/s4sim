#/bin/bash

# split0 = split3 = split6 = split7 = split8

toast_ground_schedule \
    @schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --sun-avoidance-angle 90 \
    --moon-avoidance-angle 45 \
    --sun-avoidance-altitude-deg 5 \
    --out schedule_sat.split0.txt \
    >& schedule_sat.split0.log &

# split1 = split4

toast_ground_schedule \
    @schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --sun-avoidance-angle 60 \
    --moon-avoidance-angle 45 \
    --sun-avoidance-altitude-deg 5 \
    --out schedule_sat.split1.txt \
    >& schedule_sat.split1.log &

# split2 = split5 = split9

toast_ground_schedule \
    @schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --sun-avoidance-angle 45 \
    --moon-avoidance-angle 45 \
    --sun-avoidance-altitude-deg 5 \
    --out schedule_sat.split2.txt \
    >& schedule_sat.split2.log &
