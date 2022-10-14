#/bin/bash

export TOAST_LOGLEVEL=DEBUG

# Run all the schedules

toast_ground_schedule \
    @schedule_sat.par \
    >& get_schedule.log &

#toast_ground_schedule \
#    @schedule_sat.7month.par \
#    >& get_schedule.7month.log &
