#!/bin/bash

pwv_limit=3

# for flavor in lat_delensing_sun90bk lat_wide_supplement lat_roman_supplement lat_wide; do
for flavor in lat_delensing_sun90bk; do
    echo $flavor
    schedule=schedule.txt

    case $flavor in
        lat_delensing_sun90bk)
            fname1="../scan_strategy/lat_delensing_sun90bk/schedule_lat_delensing_sun90bk.${pwv_limit}mm.break.txt"
            fname2="../scan_strategy/lat_delensing_sun90bk/schedule_lat_delensing_sun90bk.${pwv_limit}mm.season.txt"
            ;;
        lat_wide)
            fname1="../scan_strategy/lat_wide/schedule_lat_wide.${pwv_limit}mm.break.txt"
            fname2="../scan_strategy/lat_wide/schedule_lat_wide.${pwv_limit}mm.season.txt"
            ;;
        lat_wide_supplement)
            fname1="../scan_strategy/lat_wide/schedule_lat_wide_supplement.${pwv_limit}mm.break.txt"
            fname2="../scan_strategy/lat_wide/schedule_lat_wide_supplement.${pwv_limit}mm.season.txt"
            ;;
        lat_roman_supplement)
            fname1="../scan_strategy/lat_roman/schedule_lat_roman_supplement.${pwv_limit}mm.break.txt"
            fname2="../scan_strategy/lat_roman/schedule_lat_roman_supplement.${pwv_limit}mm.season.txt"
            ;;
    esac

    cat $fname1 > $schedule
    let nline=$(cat $fname2 | wc -l)
    let nline-=3
    tail -n $nline $fname2 >> $schedule

    split_dir=split_schedules/${flavor}

    # Split the observing schedule by start date of each observation

    mkdir -p ${split_dir}
    rm -f ${split_dir}/*
    let nline_tot=$(cat $schedule | wc -l)
    let nline=nline_tot-3
    echo $nline_tot $nline
    head -n 3 $schedule > header.txt
    tail -n $nline $schedule > body.txt
    for start in $(awk '{print $1}' body.txt | uniq); do
        echo $start
        schedule_out="${split_dir}/split_schedule_${start}.txt"
        cat header.txt > $schedule_out
        # grep -E "^  $start" body.txt >> $schedule_out
        grep -E "^ $start" body.txt >> $schedule_out
    done
    
    rm header.txt
    rm body.txt
done
