#!/bin/bash
#SBATCH --partition=debug
###SBATCH --qos=premium
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --job-name=CMBS4_phase1
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107a

ulimit -c unlimited
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export TOAST_FUNCTIME=1
export HDF5_USE_FILE_LOCKING=FALSE

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
let groupsize=1

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"

bands=( f020 f030 f040 f090 f150 f220 f280 )
fproot=../focalplanes/focalplane_LAT0_CHLAT

for flavor in lat_wide lat_delensing lat_delensing_core; do
    args=()
    args+="--config params.lat.toml"
    args+=" --telescope LAT"
    args+=" --job_group_size ${groupsize}"
    case $flavor in
	lat_wide)
            schedule_in="../scan_strategy/lat_wide/schedule_lat_wide.txt"
            args+=" --sim_ground.scan_cosecant_modulation"
            scan_rate="Quantity('0.5 deg / s')"
	    ;;
	lat_delensing)
            schedule_in="../scan_strategy/lat_wide/schedule_lat_delensing.txt"
            args+=" --sim_ground.scan_rate_az Quantity('1.0 deg / s')"
	    ;;
	lat_delensing_core)
            schedule_in="../scan_strategy/lat_wide/schedule_lat_delensing_core.txt"
            args+=" --sim_ground.scan_rate_az Quantity('1.0 deg / s')"
	    ;;
    esac
    schedule=$(basename $schedule_in)
    schedule=${schedule/.txt/.sorted.txt}
    if [[ ! -e $schedule ]]; then
        echo "Sorting $schedule_in to make $schedule"
        python sort_schedule.py $schedule_in $schedule >& ${schedule/.txt/.log}
    fi
    args+=" --schedule $schedule"
    for band in ${bands[*]}; do
	case $band in
	    f020|f030|f040)
                pwv_limit=3.0
		;;
	    f220|f280)
                pwv_limit=2.0
                args+=" --thinfp 8"
		;;
	    *)
                pwv_limit=2.0
                args+=" --thinfp 8"
		;;
	esac
	args+=" --focalplane ${fproot}_${band}.h5"
        args+=" --pwv_limit ${pwv_limit}"
        args+=" --sample_rate 1"
        args+=" --pixels_healpix_radec.nside 512"

	outdir=outputs/lat/$flavor/$band
        args+=" --out $outdir"
	mkdir -p $outdir
	logfile=$outdir/log.${flavor}.${band}

	if [[ -e $outdir/mapmaker_cov.fits ]]; then
	    echo "$flavor $band already complete! Skipping..."
	    continue
	fi

	echo "Writing to $logfile"
        echo "${args[@]}"
        srun -N $nnode -n $ntask -c $ncore --cpu_bind=cores \
             toast_sim_ground.py ${args[@]} \
             --sim_ground.scan_rate_az "${scan_rate}" \
             >& $logfile
        # --sim_ground.scan_rate_az "Quantity('0.5 deg / s')" \
        exit
    done
done
