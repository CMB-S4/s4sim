#!/bin/bash
#SBATCH --qos=preempt
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --job-name=CMBS4_DC0_spsat_lf_obsmat
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107

export TOAST_LOGLEVEL=DEBUG

# Perlmutter-specific fixes
export FI_CXI_OPTIMIZED_MRS="false"
export MPI4PY_RC_RECV_MPROBE="False"

# Python environment
ulimit -c unlimited
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export HDF5_USE_FILE_LOCKING=FALSE

# TOAST variables
export TOAST_FUNCTIME=1

# The number of process groups must be divisible
# by the number wafers: 12

# Parallelization
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=1
let ntask_group=$nnode_group*$ntask_node
#let groupsize=nnode_group*ntask_node
let groupsize=$ntask_node
let ngroup=$nnode/$nnode_group

echo "$(date) : Start"
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "      ntask_group = ${ntask_group}"
echo "        groupsize = ${groupsize}"
echo "           ngroup = ${ngroup}"

telescope=spsat
site=pole

for band in f030 f040; do
    case $band in
        f030|f040)
            TELESCOPE=SAT3_SAT
            tube=ST8
            ;;
        *)
            echo "Unknown band: $band"
            exit
            ;;
    esac

    echo "$(date) : Listing schedules"
    
    fnames=(`python3 ../get_fnames_serial.py split_schedule logs/${TELESCOPE} ${TELESCOPE} ${band}`)
    echo "$(date) : Found ${#fnames[@]} schedule files"
    
    # Random wait time to reduce clashes
    sleep $((RANDOM % 15))

    echo "$(date) : Looking for schedule at" `date`

    for schedule in ${fnames[*]}; do
        rootname=`basename $schedule .txt`
        logdir=logs/${TELESCOPE}/${band}
        mkdir -p $logdir
        logfile=$logdir/${rootname}.log
        logfile2=cleared_$logfile
        if [[ ! -e $logfile && ! -e $logfile2 ]]; then
            echo "$(date) : Writing $logfile"
            date > ${logfile}
            
            outdir=outputs/${TELESCOPE}/${band}/${rootname}
            mkdir -p --mode=g+rXs-w $outdir
            chgrp -R cmbs4 $outdir
            
            srun -N $nnode_group -n $ntask_group -c $ncore --cpu_bind=cores \
                 toast_sim_ground.py \
                 --config ../params/common.toml ../params/scanning_${telescope}.toml \
                 ../params/atmosphere_${site}.toml ../params/reduce_${telescope}.toml \
                 --focalplane ../focalplanes/focalplane_${TELESCOPE}_${band}.h5 \
                 --telescope $TELESCOPE \
                 --schedule $schedule \
                 --out $outdir \
                 --job_group_size ${groupsize} \
                 --mem_count.enable \
                 --scan_healpix_map.disable \
                 --pixels_healpix_radec.nside 512 \
                 --sim_noise.disable \
                 --sim_atmosphere.disable \
                 --raw_statistics.disable \
                 --filtered_statistics.disable \
                 --filterbin.report_memory \
                 --filterbin.name filterbin_${rootname} \
                 --filterbin.no_write_hits \
                 --filterbin.noiseweight_obs_matrix \
                 --filterbin.write_invcov \
                 --filterbin.write_obs_matrix \
                 --processing_mask.file ../foreground_sim/input_maps/mask_mediumcomplexity.${telescope}.${band}.h5 \
                 --save_hdf5.disable \
                 --sim_ground.median_weather \
                 --yield_cut.disable \
                 >> ${logfile} 2>&1
            echo "$(date) : Job completed"
            exit
        else
            echo "$(date) : $logfile exists"
        fi
    done
done

echo "$(date) : Nothing to do"
