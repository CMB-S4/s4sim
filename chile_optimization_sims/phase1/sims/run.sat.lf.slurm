#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --job-name=CMBS4_phase1
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
let groupsize=1

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"

flavor=sat
# for band in f030 f040 f085 f095 f145 f155 f220 f280; do
for band in f030 f040; do
    args=()
    args+="--config params.sat.toml"
    args+=" --telescope SAT"
    args+=" --job_group_size ${groupsize}"
    schedule="../scan_strategy/sat/schedule_${flavor}.txt"
    args+=" --schedule $schedule --sort_schedule"
    case $band in
	f030|f040)
	    fpfile=../focalplanes/focalplane_SAT3_SAT_${band}_ST8.h5
	    pwv_limit=3.0
            # args+=" --thinfp 16"
	    ;;
	f085|f145)
	    fpfile=../focalplanes/focalplane_SAT1_SAT_${band}_ST0.h5
	    pwv_limit=3.0
            # args+=" --thinfp 16"
	    ;;
	f095|f155)
	    fpfile=../focalplanes/focalplane_SAT1_SAT_${band}_ST1.h5
	    pwv_limit=3.0
            # args+=" --thinfp 16"
	    ;;
	f220|f280)
	    fpfile=../focalplanes/focalplane_SAT1_SAT_${band}_ST1.h5
	    pwv_limit=2.0
            # args+=" --thinfp 16"
	    ;;
        *)
            echo "Unknown band: $band"
            exit
    esac
    args+=" --focalplane ${fpfile}"
    args+=" --pwv_limit ${pwv_limit}"
    args+=" --sample_rate 1"
    
    outdir=outputs/${flavor}/$band
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
