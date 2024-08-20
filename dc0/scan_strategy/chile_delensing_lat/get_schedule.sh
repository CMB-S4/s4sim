#/bin/bash

if [[ ! -e patches_lat.txt ]]; then
    python3 make_lat_tiles.py
fi

export TOAST_LOGLEVEL=DEBUG

# Run all the schedules

toast_ground_schedule \
    @schedule_lat.par \
    @patches_lat.txt \
    --block-out 01/01-04/01 \
    --out schedules/delens.txt \
    >& get_schedule.delens.log &

toast_ground_schedule \
    @schedule_lat.par \
    --block-out 01/01-04/01 \
    --patch RISING_SCAN_40,HORIZONTAL,1.00,30.00,150.00,40.00,1440 \
    --patch SETTING_SCAN_40,HORIZONTAL,1.00,210.00,330.00,40.00,1440 \
    --out schedules/wide.txt \
    >& get_schedule.wide.log &

#for priority in 0.10 1.00 10.0 100 1000 10000 100000 1000000; do; do
for priority in 10000000; do

    toast_ground_schedule \
	@schedule_lat.par \
	@patches_lat.txt \
	--block-out 01/01-04/01 \
	--patch RISING_SCAN_40,HORIZONTAL,$priority,30.00,150.00,40.00,1440 \
	--patch SETTING_SCAN_40,HORIZONTAL,$priority,210.00,330.00,40.00,1440 \
	--out schedules/priority_$priority.txt \
	>& get_schedule.$priority.log &
    
done
