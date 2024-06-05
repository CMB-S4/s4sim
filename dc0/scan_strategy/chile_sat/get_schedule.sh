#/bin/bash

if [[ ! -e patches_sat.txt ]]; then
    python3 make_sat_tiles.py
fi

# export TOAST_LOGLEVEL=DEBUG

# Run all the schedules

logfile=get_schedule.log
echo "Writing $logfile"

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --block-out 01/01-04/01 \
    >& $logfile &
