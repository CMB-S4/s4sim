#!/bin/bash
#SBATCH --qos=preempt
#SBATCH --time=12:00:00
#SBATCH --nodes=4
#SBATCH --job-name=CMBS4_DC0_chlat_lf_noise
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107

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

# Parallelization
export OMP_NUM_THREADS=2
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
let groupsize=ntask_node
let ngroup=$nnode/$nnode_group

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "      ntask_group = ${ntask_group}"
echo "        groupsize = ${groupsize}"
echo "           ngroup = ${ngroup}"

telescope=chlat
site=chile
TELESCOPE=LAT0_CHLAT
#bands=( f030 )
bands=(f030 f040)
#bands=(f090 f150)
#bands=(f090 f150 f220 f280)

echo "Listing schedules at" `date`

fnames1=(`python3 ../get_fnames_serial.py ../split_schedules_1_upto2mm/${telescope} logs/${TELESCOPE} ${TELESCOPE} ${bands[*]}`)
fnames2=(`python3 ../get_fnames_serial.py ../split_schedules_1_over2mm/${telescope} logs/${TELESCOPE} ${TELESCOPE} ${bands[*]}`)
fnames=( ${fnames1[@]} ${fnames2[@]} )
echo "Found ${#fnames[@]} schedule files"

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

echo "Looking for schedule at" `date`

for band in ${bands[*]}; do
    logdir=logs/${TELESCOPE}/${band}
    mkdir -p $logdir
done

for band in ${bands[*]}; do
    for schedule in ${fnames[*]}; do
        rootname=`basename $schedule .txt`
        logdir=logs/${TELESCOPE}/${band}
        mkdir -p $logdir
        logfile=$logdir/${rootname}.log
        logfile2=cleared_$logdir/${rootname}.log
        logfile3=staged_$logdir/${rootname}.log
        if [[ ! -e $logfile && ! -e $logfile2 && ! -e $logfile3 ]]; then
            echo "Writing $logfile at" `date`
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
                 --sim_atmosphere.cache_dir atm_cache_${TELESCOPE} \
                 --mapmaker.report_memory \
                 --mapmaker.name mapmaker_${rootname} \
                 --mapmaker.write_invcov \
                 --mapmaker.write_hits \
                 --processing_mask.file ../foreground_sim/input_maps/mask_mediumcomplexity.${telescope}.${band}.h5 \
                 --save_hdf5.enable \
                 --save_hdf5.volume $outdir \
                 --sim_ground.session_split_key tube \
                 >> ${logfile} 2>&1 &
            sleep 1
            # Wait until the srun process has time to spawn another one
            # Running jobs always have two srun instances
            let nsrun=`ps | grep srun | wc -l`
            let running=${nsrun}/2
            let nsrun_estimate=${running}*2
            while [[ $nsrun -ne $nsrun_estimate ]]; do
                echo "There are $running running processes and $nsrun srun instances at" `date`
                echo "Waiting for 15 seconds."
                sleep 15
                let nsrun=`ps | grep srun | wc -l`
                let running=${nsrun}/2
                let nsrun_estimate=${running}*2
            done
            # Now see if we have all the slots taken
            while [[ $running -ge $ngroup ]]; do
                # We do. Wait for at least one slot to be freed
                echo "There are $running running processes at" `date`
                echo "Waiting for 15 seconds."
                sleep 15
                let running=`ps | grep srun | wc -l`/2
            done
            echo "Finding more work at" `date`
        else
            echo "$logfile exists"
        fi
    done
done

echo "Waiting for $running jobs to complete at" `date`

wait

echo "Jobs completed at" `date`
