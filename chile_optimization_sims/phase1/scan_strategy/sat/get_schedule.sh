#/bin/bash

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --out schedule_sat.txt \
    >& schedule_sat.log
