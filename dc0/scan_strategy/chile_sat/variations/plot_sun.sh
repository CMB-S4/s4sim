#/bin/bash

for month in {1..12}; do
    for angle in 30 45 60 90; do
	smonth1=$(printf "%02i" $month)
	if [[ $month == 12 ]]; then
	    year=2028
	    month2=1
	else
	    year=2027
	    let month2=$month+1
	fi
	smonth2=$(printf "%02i" $month2)
	toast_ground_schedule \
	    @chile_schedule_sat.par \
	    @patches_sat.txt \
	    --out schedules/temp.txt \
	    --debug \
	    --start "2027-$smonth1-01 00:00:00" \
	    --stop "$year-$smonth2-01 00:00:00" \
	    --moon-avoidance-angle 0 \
	    --sun-avoidance-angle $angle \
	    >& log.temp
	mv patches.png patches$smonth1.${angle}deg.png
    done
done
#	--polmap foregrounds_100GHz_fwhm1deg.fits \
#	--pol-max 1000 \
#    mv patches.png patches$smonth1.polmap.png
