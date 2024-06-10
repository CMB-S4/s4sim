#!/bin/bash

schedule_in="schedules/sidereal.txt"
schedule_out="schedules/sample.txt"
head -n 3 $schedule_in > $schedule_out
grep -E '^ 2027-06-23|^ 2027-06-24|^ 2027-06-25|^ 2027-06-26' $schedule_in \
     | grep -v "North_Direction" \
     | grep -v "Tier3" \
     | grep -v "DEC+" \
     >> $schedule_out
