#!/bin/bash
#SBATCH --qos=premium
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --job-name=CMBS4_phase3
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107a

# Python environment
ulimit -c unlimited
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export HDF5_USE_FILE_LOCKING=FALSE

# TOAST variables
export TOAST_FUNCTIME=1

# Parallelization
export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS
let groupsize=4
let ngroup=$ntask/$groupsize

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"
echo "           ngroup = ${ngroup}"

# bands=( f020 f030 f040 f090 f150 f220 f280 )
# bands=( f090 f150 f220 f280 )
bands=( f090 f150 )
fproot=../../phase2/focalplanes/focalplane_LAT0_CHLAT

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

#for flavor in lat_wide lat_delensing lat_delensing_core lat_delensing_tiled; do
#for flavor in lat_wide lat_delensing_tiled lat_delensing_core; do
# for flavor in lat_delensing_sun90bk; do
for flavor in lat_delensing_sun90bk lat_wide_supplement lat_roman_supplement lat_wide; do
# for flavor in lat_wide; do
    args=()
    args+="--config params.lat.daily.toml"
    args+=" --telescope LAT"
    args+=" --job_group_size ${groupsize}"
    for band in ${bands[*]}; do
	case $band in
	    f020|f030|f040)
                pwv_limit=3
		;;
	    f220|f280)
                pwv_limit=2
                args+=" --thinfp 16"
		;;
	    *)
                pwv_limit=3
                args+=" --thinfp 16"
		;;
	esac
        case $flavor in
	    lat_delensing_sun90bk)
                scan_rate="Quantity('1.0 deg / s')"
	        ;;
	    lat_wide)
                args+=" --sim_ground.scan_cosecant_modulation"
                scan_rate="Quantity('0.5 deg / s')"
	        ;;
	    lat_wide_supplement)
                args+=" --sim_ground.scan_cosecant_modulation"
                scan_rate="Quantity('0.5 deg / s')"
	        ;;
	    lat_roman_supplement)
                scan_rate="Quantity('1.0 deg / s')"
	        ;;
            *)
                echo "Unknown flavor = ${flavor}"
                exit
        esac
	args+=" --focalplane ${fproot}_${band}.h5"
        args+=" --pwv_limit ${pwv_limit}"
        args+=" --sample_rate 1"
        args+=" --pixels_healpix_radec.nside 512"

        split_dir=split_schedules/${flavor}
        for schedule in ${split_dir}/*txt; do
            id=`basename $schedule .txt`
            id=${id/split_schedule_/}

            outdir=daily_outputs/${flavor}/${band}/${id}
	    mkdir -p $outdir
	    logfile=$outdir/log

	    if [[ -e $outdir/mapmaker_cov.fits ]]; then
	        echo "$flavor $band already complete! Skipping..."
	        continue
	    fi

	    if [[ -e $logfile ]]; then
	        echo "$logfile exists! Skipping..."
	        continue
	    fi

	    echo "$(date) : Writing to $logfile"
            echo "${args[@]}"
            srun -n $groupsize -c $ncore --cpu_bind=cores \
                 toast_sim_ground.py ${args[@]} \
                 --schedule $schedule \
                 --sim_ground.scan_rate_az "${scan_rate}" \
                 --weights_radec.mode "I" \
                 --out $outdir \
                 >& $logfile &

            sleep 1
            # Wait until the srun process has time to spawn another one
            # Running jobs always have two srun instances
            let nsrun=`ps | grep srun | wc -l`
            let running=${nsrun}/2
            let nsrun_estimate=${running}*2
            while [[ $nsrun -ne $nsrun_estimate ]]; do
                echo -n "$(date) : There are $running running processes "
                echo -n "and $nsrun srun instances. "
                echo "Waiting for 15 seconds."
                sleep 15
                let nsrun=$(ps | grep srun | wc -l)
                let running=${nsrun}/2
                let nsrun_estimate=${running}*2
            done
            # Now see if we have all the slots taken
            while [[ $running -ge $ngroup ]]; do
                # We do. Wait for at least one slot to be freed
                echo -n "$(date) : There are $running running processes."
                echo "Waiting for 15 seconds."
                sleep 15
                let running=$(ps | grep srun | wc -l)/2
            done
            echo "$(date) : Finding more work"

        done
    done
done

wait

echo "$(date) : All done!"
