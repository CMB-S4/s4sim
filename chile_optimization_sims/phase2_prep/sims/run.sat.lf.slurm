#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --job-name=CMBS4_phase2_prep
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
let groupsize=8

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

flavor=sat
for isplit in {0..9}; do
    split="split${isplit}"
    case $isplit in
        0|1|2)
            bands=(f030 f040)
            fpdir=../focalplanes_split0
            ;;
        3|4|5)
            bands=(f030 f040)
            fpdir=../focalplanes_split3
            ;;
        6)
            bands=(f030 f040)
            fpdir=../focalplanes_split6
            ;;
        7)
            bands=(f030 f040)
            fpdir=../focalplanes_split7
            ;;
        8|9)
            bands=(f030 f040)
            fpdir=../focalplanes_split8
            ;;
        *)
            echo "Unknown split: $isplit"
            exit
    esac
    case $isplit in
        0|3|6|7|8)
            schedule=../scan_strategy/sat/schedule_sat.split0.txt
            ;;
        1|4)
            schedule=../scan_strategy/sat/schedule_sat.split1.txt
            ;;
        2|5|9)
            schedule=../scan_strategy/sat/schedule_sat.split2.txt
            ;;
        *)
            echo "Unknown split: $isplit"
            exit
    esac
    for band in ${bands[*]}; do
        args=()
        args+="--config params.sat.toml"
        args+=" --telescope SAT"
        args+=" --job_group_size ${groupsize}"
        args+=" --schedule $schedule --sort_schedule"
        case $band in
	    f030|f040)
	        fpfile_in=${fpdir}/focalplane_SAT3_SAT_${band}_ST8.h5
	        pwv_limit=3.0
	        ;;
	    f085|f090|f145|f150)
	        fpfile_in=${fpdir}/focalplane_SAT1_SAT_${band}_ST0.h5
	        pwv_limit=3.0
	        ;;
	    f095|f155)
	        fpfile_in=${fpdir}/focalplane_SAT1_SAT_${band}_ST1.h5
	        pwv_limit=3.0
	        ;;
	    f220|f280)
	        fpfile_in=${fpdir}/focalplane_SAT1_SAT_${band}_ST2.h5
	        pwv_limit=2.0
	        ;;
            *)
                echo "Unknown band: $band"
                exit
        esac
        # Avoid clashes in focalplane filenames
        fpfile=focalplane_SAT_${split}_${band}.h5
        rsync $fpfile_in $fpfile
        args+=" --focalplane ${fpfile}"
        args+=" --pwv_limit ${pwv_limit}"
        args+=" --sample_rate 1"
        
        outdir=outputs/${flavor}/${split}/$band
        args+=" --out $outdir"
        mkdir -p $outdir
        logfile=$outdir/log
        
        if [[ -e $outdir/mapmaker_cov.fits ]]; then
	    echo "$flavor $split $band already complete! Skipping..."
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
             --binner.sync_type "allreduce" \
             >& $logfile
        
        exit
    done
done
