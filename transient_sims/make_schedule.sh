#!/bin/bash

# Run the TOAST scheduler

schedule="schedule.txt"
echo "Writing $schedule"
toast_ground_schedule \
    @chile_schedule_lat.par \
    --start "2030-06-15 00:00:00" \
    --stop "2030-06-16 00:00:00" \
    --out $schedule

# Break the observing schedule into separate observations

nline=1
echo "Splitting $schedule into files with $nline entries"
outdir=split_schedule
echo "Writing to $outdir"
mkdir -p ${outdir}
rm -rf ${outdir}/*

schedule_out="${outdir}/split_schedule_"
let nces=`wc -l ${schedule} | awk '{print $1 - 3}'`
echo "NCES = ${nces}"
[[ $nces -eq 0 ]] && continue

head -n 3 $schedule > header.txt
awk '{if (NR > 3) print}' $schedule | \
    split --lines=${nline} --numeric-suffixes \
          --suffix-length=4 --additional-suffix=.txt - \
          $schedule_out
for fname_in in ${schedule_out}*; do
    mv ${fname_in} temp.txt
    fname_out=`awk '{print $6 "-" $10 "-" $11}' temp.txt`
    if [[ $fname_out == CALIBRATION_BREAK* ]]; then
        # We don't simulate calibration breaks
        continue
    fi
    fname_out=$(dirname $fname_in)/${fname_out}.txt
    cat header.txt temp.txt > ${fname_out}
    rm temp.txt
done

rm header.txt
