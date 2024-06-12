#!/bin/bash

export OMP_NUM_THREADS=1
export TOAST_FUNCTIME=1

schedule=schedule.txt
head -n 4 sample.txt > $schedule

for flavor in hwpss; do
    logfile=${flavor}.log
    if [[ -e $logfile ]]; then
	echo "$logfile exists. Skipping..."
	#continue
    fi

    echo "Writing to $logfile"

    outdir=outputs/$flavor
    mkdir -p $outdir

    case $flavor in
	hwpss)
	    args="--sim_hwpss.enable"
	    args+=" --sim_ground.hwp_angle 'hwp_angle'"
	    args+=" --sim_ground.hwp_rpm 120"
	    args+=" --demodulate.enable"
	    args+=" --hwpfilter.enable"
	    ;;
	*)
	    echo "Unknown flavor: $flavor"
	    continue
	    ;;
    esac

    mpirun -np 16 python toast_sim_ground.py \
	   --focalplane focalplane_SAT1_SAT_f155_ST0.h5 \
	   --telescope SAT \
	   --schedule $schedule \
	   --sample_rate 10.0 \
	   --pixels_healpix_radec.nside 512 \
	   --thinfp 1 \
	   --out $outdir \
	   --sim_atmosphere.disable \
	   --sim_noise.disable \
	   --mapmaker.no_write_map \
	   --mapmaker.no_write_rcond \
	   --baselines.disable \
	   --polyfilter1D.disable \
	   ${args} \
	   >& $logfile

done

#	   --job_group_size 1
#	   --mapmaker.no_write_binmap
#	   --pwv_limit 3.0
#	   --weights_azel.mode "I"
#	   --weights_radec.mode "I"
