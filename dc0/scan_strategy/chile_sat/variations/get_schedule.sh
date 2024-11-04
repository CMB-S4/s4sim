#/bin/bash

# Run all the schedules

# Twilight types:
# -18 deg : Astronomical twilight (scheduler default)
# -12 deg : nautical twilight
# -6 deg : civil twilight
# 0 deg : sunrise

for avoid_el in -18 -12 -6 0 1 2 3 4 5 10 15; do
    toast_ground_schedule \
	@chile_schedule_sat.par \
	@patches_sat.txt \
	--block-out 01/01-04/01 \
	--out schedules/solar90.${avoid_el}.txt \
	--sun-avoidance-angle 90 \
	--sun-avoidance-altitude-deg " ${avoid_el}" \
	>& log.solar90.$avoid_el &
done

exit

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --block-out 01/01-04/01 \
    --out schedules/baseline.txt \
    >& log.baseline

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --out schedules/no_break.txt \
    >& log.no_break

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --block-out 01/01-04/01 \
    --moon-avoidance-angle 0 \
    --out schedules/no_lunar_avoidance.txt \
    >& log.no_lunar_avoidance

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.south_only.txt \
    --block-out 01/01-04/01 \
    --out schedules/south_only.txt \
    >& log.south_only

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --block-out 01/01-04/01 \
    --el-min 40 \
    --out schedules/el_min_40.txt \
    >& log.el_min_40

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --patch North_Direction,SIDEREAL,100.0,330.00,30.00,67.00,130,160,15 \
    --block-out 01/01-04/01 \
    --out schedules/sidereal.txt \
    >& log.sidereal

toast_ground_schedule \
    @chile_schedule_sat.par \
    @patches_sat.txt \
    --patch South_Direction,SIDEREAL,10.0,150.00,210.00,67.00,30,60,15 \
    --patch North_Direction,SIDEREAL,100.0,330.00,30.00,67.00,130,160,15 \
    --moon-avoidance-angle 0 \
    --out schedules/all.txt \
    --el-min 40 \
    >& log.all
