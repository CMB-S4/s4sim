#/bin/bash

#export TOAST_LOGLEVEL=DEBUG

for flavor in sun90max sun90ra0 sun90bk; do
    python make_targets.ra0.py ../sat/hits_season_${flavor}.fits patches.${flavor}.txt

    toast_ground_schedule \
	@schedule_lat_delensing.par \
	--start "2030-01-01 00:00:00" \
	--stop "2031-01-01 00:00:00" \
	@patches.${flavor}.txt \
	--time-step 1200 \
	--out schedule_lat_delensing_${flavor}.txt \
	>& schedule_lat_delensing_${flavor}.log &
done

# python gapfill_schedule.py schedule_sat.sun90max.txt schedule_sat.sun45max.txt schedule_sat.sun45_supplement.txt
