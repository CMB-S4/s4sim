#!/bin/bash

export OMP_NUM_THREADS=1
export TOAST_FUNCTIME=1

for flavor in sample; do
    outdir=outputs/$flavor
    mkdir -p $outdir

    logfile=$outdir/${flavor}.log
    if [[ -e $logfile ]]; then
	echo "$logfile exists. Skipping..."
	#continue
    fi

    echo "Writing to $logfile"

    mpirun -np 16 toast_sim_ground.py \
	   --focalplane ../../../focalplanes/focalplane_SAT1_SAT_f155.h5 \
	   `# --pwv_limit 3.0` \
	   --telescope SAT \
	   --schedule schedules/${flavor}.txt \
	   --sample_rate 1.0 \
	   --pixels_healpix_radec.nside 512 \
	   --thinfp 1 \
	   --out $outdir \
	   --job_group_size 1 \
	   --sim_atmosphere.disable \
	   --sim_noise.disable \
	   --mapmaker.no_write_binmap \
	   --mapmaker.no_write_map \
	   --mapmaker.no_write_rcond \
	   --baselines.disable \
	   --weights_azel.mode "I" \
	   --weights_radec.mode "I" \
	   --polyfilter1D.disable \
	   >& $logfile
done
