#/bin/bash

#export TOAST_LOGLEVEL=DEBUG

flavors=(sun90max sun90ra0 sun90bk)

#for flavor in ${flavors[*]}; do
#    python make_targets.south.py ../sat/hits_season_${flavor}.fits patches.${flavor}.txt

#    toast_ground_schedule \
#	@schedule_lat_delensing.par \
#	--start "2030-01-01 00:00:00" \
#	--stop "2031-01-01 00:00:00" \
#	@patches.${flavor}.txt \
#	--time-step 1200 \
#	--out schedule_lat_delensing_${flavor}.txt \
#	>& schedule_lat_delensing_${flavor}.log &
#done

wait

for flavor in ${flavors[*]}; do
    python gapfill_schedule.py \
	   schedule_lat_delensing_${flavor}.txt \
	   ../lat_wide/schedule_lat_wide.txt \
	   schedule_lat_delensing_${flavor}+wide.txt \
	   merge

    python gapfill_schedule.py \
	   schedule_lat_delensing_${flavor}.txt \
	   ../lat_wide/schedule_lat_wide.txt \
	   schedule_lat_delensing_${flavor}_gapfill.txt
done
