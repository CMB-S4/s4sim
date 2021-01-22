#!/bin/bash

for site in pole chile; do
    echo $site
    for tel in lat sat; do
        echo $tel
        nline=1

        indir=scan_strategy/${site}_${tel}/schedules
        outdir=scan_strategy/${site}_${tel}/split_schedules
        mkdir -p ${outdir}
        rm -rf ${outdir}/*

        schedule_in="${indir}/${site}_schedule_${tel}.txt"
        schedule_out="${outdir}/split_schedule_${tel}_"
        let nces=`wc -l ${schedule_in} | awk '{print $1 - 3}'`
        echo "NCES = ${nces}"
        head -n 3 $schedule_in > header.txt
        awk '{if (NR > 3) print}' $schedule_in | \
            split --lines=${nline} --numeric-suffixes \
                  --suffix-length=4 --additional-suffix=.txt - \
                  $schedule_out
        for fname in ${schedule_out}*; do
            mv ${fname} temp.txt
            cat header.txt temp.txt > ${fname}
            rm temp.txt
        done
    done
done
