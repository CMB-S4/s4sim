#/bin/bash

export TOAST_LOGLEVEL=DEBUG

# Run all the schedules

toast_ground_schedule \
    @schedule_lat.par \
    >& get_schedule.log &
