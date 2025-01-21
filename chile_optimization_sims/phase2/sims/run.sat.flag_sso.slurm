#!/bin/bash
#SBATCH --qos=premium
#SBATCH --time=00:30:00
#SBATCH --nodes=8
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

flavors=(sun45)
for flavor in ${flavors[*]}; do
    for band in f155; do
        for period in season break; do
            args=()
            args+="--config params.sat.toml"
            args+=" --telescope SAT"
            args+=" --job_group_size ${groupsize}"
            case $band in
	        f030|f040)
	            fpfile=../focalplanes/focalplane_SAT3_SAT_${band}_ST8.h5
	            pwv_limit=3
                    # args+=" --thinfp 16"
	            ;;
	        f085|f145|f090|f150)
	            fpfile=../focalplanes/focalplane_SAT1_SAT_${band}_ST0.h5
	            pwv_limit=3
                    # args+=" --thinfp 16"
	            ;;
	        f095|f155)
	            fpfile=../focalplanes/focalplane_SAT1_SAT_${band}_ST1.h5
	            pwv_limit=3
                    # args+=" --thinfp 16"
	            ;;
	        f220|f280)
	            fpfile=../focalplanes/focalplane_SAT1_SAT_${band}_ST2.h5
	            pwv_limit=2
                    # args+=" --thinfp 16"
	            ;;
                *)
                    echo "Unknown band: $band"
                    exit
            esac
            schedule="../scan_strategy/sat/schedule_sat.${flavor}.${pwv_limit}mm.${period}.txt"
            args+=" --schedule $schedule --sort_schedule"
            args+=" --focalplane ${fpfile}"
            args+=" --pwv_limit ${pwv_limit}"
            args+=" --sample_rate 1"

            for radius in 90 75 60 45; do
                full_flavor=${flavor}_flag$radius
                outdir=outputs/${full_flavor}/$band/${period}
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
                     --binner.sync_type "allreduce" \
                     --flag_sso.enable \
                     --flag_sso.sso_names "['Sun']" \
                     --flag_sso.sso_radii "[Quantity('${radius} deg')]" \
                     >& $logfile

                exit
            done
        done
    done
done
