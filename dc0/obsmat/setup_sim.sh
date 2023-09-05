#!/bin/bash

python3 find_irreducible_scanset.py ../scan_strategy/pole_sat/schedules/pole_schedule_sat.pruned.upto2mm.txt --pole_mode --tol 0.05

# Split schedules for batch processing

nline=1
telescope=spsat
schedule_in=irreducible_schedule.txt
if [[ ! -e $schedule_in ]]; then
    echo "No such schedule: $schedule_in"
    exit
else
    echo "Splitting $schedule_in"
fi

outdir=split_schedule
echo "Writing to $outdir"
mkdir -p ${outdir}
rm -rf ${outdir}/*

schedule_out="${outdir}/split_schedule_"
let nces=`wc -l ${schedule_in} | awk '{print $1 - 3}'`
echo "NCES = ${nces}"
[[ $nces -eq 0 ]] && continue

head -n 3 $schedule_in > header.txt
awk '{if (NR > 3) print}' $schedule_in | \
    split --lines=${nline} --numeric-suffixes \
          --suffix-length=4 --additional-suffix=.txt - \
          $schedule_out

for fname_in in ${schedule_out}*; do
    mv ${fname_in} temp.txt
    fname_out=`awk '{print $6 "-" $9 "-" $5}' temp.txt`
    fname_out=$(dirname $fname_in)/${fname_out}.txt
    cat header.txt temp.txt > ${fname_out}
    rm temp.txt
done

rm header.txt
