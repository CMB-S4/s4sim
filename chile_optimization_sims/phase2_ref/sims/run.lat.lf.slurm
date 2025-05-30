#!/bin/bash
#SBATCH --qos=premium
#SBATCH --time=02:00:00
#SBATCH --nodes=1
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
bands=( f020 f030 f040 )
fproot=../focalplanes/focalplane_LAT2_SPLAT

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

#for flavor in lat_wide lat_delensing lat_delensing_core lat_delensing_tiled; do
for flavor in splat; do
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
	        splat)
                    schedule="../scan_strategy/splat/schedule_splat.${pwv_limit}mm.${period}.txt"
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
                 >& $logfile
            exit
        done
    done
done
