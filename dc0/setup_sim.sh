#!/bin/bash

# Generate TOAST focalplane files

if [[ 1 -eq 0 ]]; then

mkdir -p focalplanes
cd focalplanes
s4_hardware_to_toast3.py --telescope LAT0
s4_hardware_to_toast3.py --telescope LAT1
s4_hardware_to_toast3.py --telescope LAT2
s4_hardware_to_toast3.py --telescope SAT1
s4_hardware_to_toast3.py --telescope SAT2
s4_hardware_to_toast3.py --telescope SAT3

fi

# Split schedules for batch processing

for suffix in .upto2mm .over2mm; do
#for suffix in .upto2mm_with_break .upto3mm_with_break; do
    for nline in 1; do
        # for telescope in chlat splat spsat; do
        for telescope in splat spsat; do
            case $telescope in
                chlat)
                    schedule_in=scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned${suffix}.txt
                    ;;
                splat)
                    schedule_in=scan_strategy/pole_lat/schedules/pole_schedule_lat.pruned${suffix}.txt
                    ;;
                spsat)
                    schedule_in=scan_strategy/pole_sat/schedules/pole_schedule_sat.pruned${suffix}.txt
                    ;;
            esac

            if [[ ! -e $schedule_in ]]; then
                echo "No such schedule: $schedule_in"
                continue
            else
                echo "Splitting $schedule_in"
            fi

            outdir=split_schedules_${nline}${suffix/./_}/${telescope}
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
                if [[ $telescope == "chlat" ]]; then
                    fname_out=`awk '{print $8 "-" $22 "-" $23}' temp.txt`
                else
                    fname_out=`awk '{print $6 "-" $10 "-" $11}' temp.txt`
                fi
                if [[ $fname_out == CALIBRATION_BREAK* ]]; then
                    # We don't simulate calibration breaks
                    continue
                fi
                fname_out=$(dirname $fname_in)/${fname_out}.txt
                cat header.txt temp.txt > ${fname_out}
                rm temp.txt
            done

            rm header.txt
        done
    done
done
