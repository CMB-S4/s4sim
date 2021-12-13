#!/bin/bash

# create fake focalplanes for caching atmosphere

toast_fake_focalplane --fov_deg 9.0 --out dummy_focalplane.440Hz.h5 --min_pix 1 --sample_rate 440
toast_fake_focalplane --fov_deg 9.0 --out dummy_focalplane.220Hz.h5 --min_pix 1 --sample_rate 220
toast_fake_focalplane --fov_deg 30.0 --out dummy_focalplane.20Hz.h5 --min_pix 1 --sample_rate 20

# split the observing schedules

for nline in 1 32; do
    for telescope in chlat splat spsat; do
        case $telescope in
            chlat)
                schedule_in=../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned.txt
                ;;
            splat)
                schedule_in=../scan_strategy/pole_lat/schedules/pole_schedule_lat.pruned.txt
                ;;
            spsat)
                schedule_in=../scan_strategy/pole_sat/schedules/pole_schedule_sat.pruned.txt
                ;;
        esac
        
        outdir=split_schedules_${nline}/${telescope}
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
        for fname in ${schedule_out}*; do
            mv ${fname} temp.txt
            cat header.txt temp.txt > ${fname}
            rm temp.txt
        done
        
        rm header.txt
    done
done