#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --job-name=CMBS4_phase2
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
export OMP_NUM_THREADS=4
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

# bands=( f020 f030 f040 f090 f150 f220 f280 )
bands=( f090 f150 f220 f280 )
fproot=../focalplanes/focalplane_LAT0_CHLAT

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

for flavor in lat_wide lat_delensing; do
    args=()
    args+="--config params.lat.toml"
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
                pwv_limit=2
                args+=" --thinfp 16"
		;;
	esac
        for period in season break; do
            case $flavor in
	        lat_wide)
                    schedule="../scan_strategy/lat_wide/schedule_lat_wide.${pwv_limit}mm.${period}.txt"
                    args+=" --sim_ground.scan_cosecant_modulation"
                    scan_rate="Quantity('0.5 deg / s')"
	            ;;
	        lat_delensing)
                    schedule="../scan_strategy/lat_delensing/schedule_lat_delensing.${pwv_limit}mm.${period}.txt"
                    scan_rate="Quantity('1.0 deg / s')"
	            ;;
            esac
            args+=" --schedule $schedule --sort_schedule"
	    args+=" --focalplane ${fproot}_${band}.h5"
            args+=" --pwv_limit ${pwv_limit}"
            args+=" --sample_rate 1"
            args+=" --pixels_healpix_radec.nside 512"
            
            outdir=outputs/${flavor}/${band}/${period}
            args+=" --out $outdir"
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
done
