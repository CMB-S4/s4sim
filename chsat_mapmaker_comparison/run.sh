#!/bin/bash

export OMP_NUM_THREADS=1
export TOAST_FUNCTIME=1

schedule=schedule.txt
head -n 4 sample.txt > $schedule

# for flavor in hwpss ps atmosphere; do
for flavor in atmosphere; do
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
	    #args+=" --demodulate.enable"
	    args+=" --hwpfilter.enable"
	    args+=" --hwpfilter.filter_order 10"
	    args+=" --hwpfilter.trend_order 5"
	    ;;
	ps)
	    args="--scan_healpix_map.enable"
	    args+=" --scan_healpix_map.file sim_ps_map_150_24arcmin_celestial.fits"
	    args+=" --pixels_healpix_radec.nside 2048"
	    args+=" --pixels_healpix_radec.no_nest"
	    ;;
	atmosphere)
	    args="--sim_atmosphere.enable"
	    args+=" --sim_atmosphere.cache_dir 'atm_cache'"
	    args+=" --sim_atmosphere.field_of_view Quantity('35deg')"
	    args+=" --sim_atmosphere.lmin_center Quantity('0.001m')"
	    args+=" --sim_atmosphere.lmin_sigma Quantity('0.0001m')"
	    args+=" --sim_atmosphere.lmax_center Quantity('1.0m')"
	    args+=" --sim_atmosphere.lmax_sigma Quantity('0.1m')"
	    args+=" --sim_atmosphere.xstep Quantity('4m')"
	    args+=" --sim_atmosphere.ystep Quantity('4m')"
	    args+=" --sim_atmosphere.zstep Quantity('4m')"
	    args+=" --sim_atmosphere.zmax Quantity('200m')"
	    args+=" --sim_atmosphere.gain 4e-5"
	    args+=" --sim_atmosphere.wind_dist Quantity('1000m')"
	    ;;
	*)
	    echo "Unknown flavor: $flavor"
	    continue
	    ;;
    esac

    OMP_NUM_THREADS=16 mpirun -np 1 \
        python toast_sim_ground.py \
	--focalplane focalplane_SAT1_SAT_f155_ST0.h5 \
	--telescope SAT \
	--schedule $schedule \
	--sample_rate 100.0 \
	--pixels_healpix_radec_final.enable \
	--pixels_healpix_radec_final.nside 512 \
	--thinfp 16 \
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
