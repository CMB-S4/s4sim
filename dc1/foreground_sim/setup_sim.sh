#!/bin/bash

for suffix in .upto2mm .over2mm; do
    for nline in 1; do
        for telescope in chlat splat spsat; do
            case $telescope in
                chlat)
                    schedule_in=../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned${suffix}.txt
                    ;;
                splat)
                    schedule_in=../scan_strategy/pole_lat/schedules/pole_schedule_lat.pruned${suffix}.txt
                    ;;
                spsat)
                    schedule_in=../scan_strategy/pole_sat/schedules/pole_schedule_sat.pruned${suffix}.txt
                    ;;
            esac
            
            outdir=split_schedules_${nline}${suffix/./_}/${telescope}
            mkdir -p ${outdir}
            rm -rf ${outdir}/*
            
            schedule_out="${outdir}/split_schedule_"
            let nces=`wc -l ${schedule_in} | awk '{print $1 - 3}'`
            echo "NCES = ${nces}"
            head -n 3 $schedule_in > header.txt
            awk '{if (NR > 3) print}' $schedule_in | \
                split --lines=${nline} --numeric-suffixes \
                      --suffix-length=4 --additional-suffix=.txt - \
                      $schedule_out
            for fname_in in ${schedule_out}*; do
                mv ${fname_in} temp.txt
                fname_out=`awk '{print $8 "-" $22 "-" $23}' temp.txt`
                fname_out=$(dirname $fname_in)/${fname_out}.txt
                cat header.txt temp.txt > ${fname_out}
                rm temp.txt
            done

            rm header.txt
        done
    done
done
