#/bin/bash

#export TOAST_LOGLEVEL=DEBUG

scantime=60

for angle in 45 90; do

toast_ground_schedule \
    @schedule_sat.par \
    --start "2030-01-01 00:00:00" \
    --stop "2031-01-01 00:00:00" \
    --sun-avoidance-angle $angle \
    `#                               weight, lon,  lat,radius,throw,scantime_min` \
    `# Southern patch` \
    --patch South_Deep11,MAX-DEPTH,1e12,3.0,-45.0,15.000,40.00,$scantime \
    --patch South_Deep12,MAX-DEPTH,1e9,3.0,-50.0,15.000,40.00,$scantime \
    --patch South_Deep13,MAX-DEPTH,1e6,3.0,-55.0,15.000,40.00,$scantime \
    --patch South_Deep14,MAX-DEPTH,1e9,3.0,-60.0,15.000,40.00,$scantime \
    --patch South_Deep15,MAX-DEPTH,1e12,3.0,-65.0,15.000,40.00,$scantime \
    --patch South_Deep16,MAX-DEPTH,1e15,3.0,-70.0,15.000,40.00,$scantime \
    `#` \
    --patch South_Deep21,MAX-DEPTH,1e9,-2.0,-45.0,15.000,40.00,$scantime \
    --patch South_Deep22,MAX-DEPTH,1e6,-2.0,-50.0,15.000,40.00,$scantime \
    --patch South_Deep23,MAX-DEPTH,1e3,-2.0,-55.0,15.000,40.00,$scantime \
    --patch South_Deep24,MAX-DEPTH,1e6,-2.0,-60.0,15.000,40.00,$scantime \
    --patch South_Deep25,MAX-DEPTH,1e9,-2.0,-65.0,15.000,40.00,$scantime \
    --patch South_Deep26,MAX-DEPTH,1e12,-2.0,-70.0,15.000,40.00,$scantime \
    `#` \
    --patch South_Deep31,MAX-DEPTH,1e6,-7.0,-45.0,15.000,40.00,$scantime \
    --patch South_Deep32,MAX-DEPTH,1e3,-7.0,-50.0,15.000,40.00,$scantime \
    --patch South_Deep33,MAX-DEPTH,1e0,-7.0,-55.0,15.000,40.00,$scantime \
    --patch South_Deep34,MAX-DEPTH,1e3,-7.0,-60.0,15.000,40.00,$scantime \
    --patch South_Deep35,MAX-DEPTH,1e6,-7.0,-65.0,15.000,40.00,$scantime \
    --patch South_Deep36,MAX-DEPTH,1e9,-7.0,-70.0,15.000,40.00,$scantime \
    `#` \
    --patch South_Deep41,MAX-DEPTH,1e9,-12.0,-45.0,15.000,40.00,$scantime \
    --patch South_Deep42,MAX-DEPTH,1e6,-12.0,-50.0,15.000,40.00,$scantime \
    --patch South_Deep43,MAX-DEPTH,1e3,-12.0,-55.0,15.000,40.00,$scantime \
    --patch South_Deep44,MAX-DEPTH,1e6,-12.0,-60.0,15.000,40.00,$scantime \
    --patch South_Deep45,MAX-DEPTH,1e9,-12.0,-65.0,15.000,40.00,$scantime \
    --patch South_Deep46,MAX-DEPTH,1e12,-12.0,-70.0,15.000,40.00,$scantime \
    `#` \
    --patch South_Deep51,MAX-DEPTH,1e12,-17.0,-45.0,15.000,40.00,$scantime \
    --patch South_Deep52,MAX-DEPTH,1e9,-17.0,-50.0,15.000,40.00,$scantime \
    --patch South_Deep53,MAX-DEPTH,1e6,-17.0,-55.0,15.000,40.00,$scantime \
    --patch South_Deep54,MAX-DEPTH,1e9,-17.0,-60.0,15.000,40.00,$scantime \
    --patch South_Deep55,MAX-DEPTH,1e12,-17.0,-65.0,15.000,40.00,$scantime \
    --patch South_Deep56,MAX-DEPTH,1e15,-17.0,-70.0,15.000,40.00,$scantime \
    `# Northern patch` \
    --patch North_Deep11,MAX-DEPTH,1e33,140.0,15.0,15.000,40.00,$scantime \
    --patch North_Deep12,MAX-DEPTH,1e30,140.0,10.0,15.000,40.00,$scantime \
    --patch North_Deep13,MAX-DEPTH,1e27,140.0,5.0,15.000,40.00,$scantime \
    --patch North_Deep14,MAX-DEPTH,1e24,140.0,0.0,15.000,40.00,$scantime \
    --patch North_Deep15,MAX-DEPTH,1e27,140.0,-5.0,15.000,40.00,$scantime \
    --patch North_Deep16,MAX-DEPTH,1e30,140.0,-10.0,15.000,40.00,$scantime \
    `#` \
    --patch North_Deep21,MAX-DEPTH,1e30,145.0,15.0,15.000,40.00,$scantime \
    --patch North_Deep22,MAX-DEPTH,1e27,145.0,10.0,15.000,40.00,$scantime \
    --patch North_Deep23,MAX-DEPTH,1e24,145.0,5.0,15.000,40.00,$scantime \
    --patch North_Deep24,MAX-DEPTH,1e21,145.0,0.0,15.000,40.00,$scantime \
    --patch North_Deep25,MAX-DEPTH,1e24,145.0,-5.0,15.000,40.00,$scantime \
    --patch North_Deep26,MAX-DEPTH,1e27,145.0,-10.0,15.000,40.00,$scantime \
    `#` \
    --patch North_Deep31,MAX-DEPTH,1e27,150.0,15.0,15.000,40.00,$scantime \
    --patch North_Deep32,MAX-DEPTH,1e24,150.0,10.0,15.000,40.00,$scantime \
    --patch North_Deep33,MAX-DEPTH,1e21,150.0,5.0,15.000,40.00,$scantime \
    --patch North_Deep34,MAX-DEPTH,1e18,150.0,0.0,15.000,40.00,$scantime \
    --patch North_Deep35,MAX-DEPTH,1e21,150.0,-5.0,15.000,40.00,$scantime \
    --patch North_Deep36,MAX-DEPTH,1e24,150.0,-10.0,15.000,40.00,$scantime \
    `#` \
    --patch North_Deep41,MAX-DEPTH,1e30,155.0,15.0,15.000,40.00,$scantime \
    --patch North_Deep42,MAX-DEPTH,1e27,155.0,10.0,15.000,40.00,$scantime \
    --patch North_Deep43,MAX-DEPTH,1e24,155.0,5.0,15.000,40.00,$scantime \
    --patch North_Deep44,MAX-DEPTH,1e21,155.0,0.0,15.000,40.00,$scantime \
    --patch North_Deep45,MAX-DEPTH,1e24,155.0,-5.0,15.000,40.00,$scantime \
    --patch North_Deep46,MAX-DEPTH,1e27,155.0,-10.0,15.000,40.00,$scantime \
    `#` \
    --patch North_Deep51,MAX-DEPTH,1e33,160.0,15.0,15.000,40.00,$scantime \
    --patch North_Deep52,MAX-DEPTH,1e30,160.0,10.0,15.000,40.00,$scantime \
    --patch North_Deep53,MAX-DEPTH,1e27,160.0,5.0,15.000,40.00,$scantime \
    --patch North_Deep54,MAX-DEPTH,1e24,160.0,0.0,15.000,40.00,$scantime \
    --patch North_Deep55,MAX-DEPTH,1e27,160.0,-5.0,15.000,40.00,$scantime \
    --patch North_Deep56,MAX-DEPTH,1e30,160.0,-10.0,15.000,40.00,$scantime \
    --time-step 900 \
    --out schedule_sat.sun${angle}bk.txt \
    >& schedule_sat.sun${angle}bk.log &
done

wait

python gapfill_schedule.py schedule_sat.sun90bk.txt schedule_sat.sun45bk.txt schedule_sat.sun45_supplement.bk.txt
