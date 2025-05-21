#/bin/bash

# This schedule targets the Northern Roman wide field
# https://roman.gsfc.nasa.gov/science/High_Latitude_Wide_Area_Survey.html

toast_ground_schedule \
    @schedule_lat.par \
    --patch roman_north_full,1,200.0,2.0,120.0,-20.000 \
    --patch roman_north_000,0.001,200.0,2.0,190.0,-20.000 \
    --patch roman_north_001,0.001,190.0,2.0,180.0,-20.000 \
    --patch roman_north_002,0.001,180.0,2.0,170.0,-20.000 \
    --patch roman_north_003,0.001,170.0,2.0,160.0,-20.000 \
    --patch roman_north_004,0.001,160.0,2.0,150.0,-20.000 \
    --patch roman_north_005,0.001,150.0,2.0,140.0,-20.000 \
    --patch roman_north_006,0.001,140.0,2.0,130.0,-20.000 \
    --patch roman_north_007,0.001,130.0,2.0,120.0,-20.000 \
    --ra-period 5 --ra-amplitude 5 \
    --out schedule_lat_roman.txt \
    >& schedule_lat_roman.log &

python gapfill_schedule.py ../lat_delensing_sun90bk/schedule_lat_delensing_sun90bk.txt schedule_lat_roman.txt schedule_lat_roman_supplement.txt
