#!/bin/bash
#SBATCH --partition=regular
#SBATCH --time=00:30:00
#SBATCH --nodes=512
#SBATCH --job-name=DC1_noise
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=knl
#SBATCH --core-spec=4
#SBATCH --account=mp107

ulimit -c unlimited
export MALLOC_MMAP_THRESHOLD_=131072
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export TOAST_FUNCTIME=1
export HDF5_USE_FILE_LOCKING=FALSE

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=1
let groupsize=nnode_group*ntask_node

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "        groupsize = ${groupsize}"

telescope=chlat
site=chile
TELESCOPE=LAT0_CHLAT
# cannot have more than 384 MPI tasks in a group to process 30 and 40GHz
bands=( f030 )
#bands=(f030 f040)
#bands=(f090 f150 f220 f280)

for band in ${bands[*]}; do
    logdir=logs/${telescope}/${band}
    mkdir -p $logdir
done

echo "Listing schedules at" `date`

#fnames=split_schedules_1/${telescope}/split_schedule_????.txt
fnames=(`python3 get_fnames.py split_schedules_1/${telescope} logs/${TELESCOPE} ${TELESCOPE} ${bands[*]}`)
echo "Found ${#fnames[@]} schedule files"

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

echo "Looking for schedule at" `date`

let running=0

for schedule in ${fnames[*]}; do
    rootname=`basename $schedule .txt`
    rootname=${TELESCOPE}_${rootname}
    for band in ${bands[*]}; do
        logdir=logs/${TELESCOPE}/${band}
        logfile=$logdir/${rootname}_${band}.log
        if [[ ! -e $logfile ]]; then
            echo "Writing $logfile at" `date`
            date > ${logfile}

            outdir=outputs/${TELESCOPE}/${band}/${rootname}
            mkdir -p $outdir
            stripe_medium $outdir
        
            srun -N 1 -n $ntask_node -c $ncore --cpu_bind=cores toast_sim_ground.py \
                 --config common.toml scanning_${telescope}.toml atmosphere_${site}.toml reduce_${telescope}.toml \
                 --focalplane ../focalplanes/focalplane_${TELESCOPE}_${band}.h5 \
                 --telescope $TELESCOPE \
                 --schedule $schedule \
                 --out $outdir \
                 --job_group_size ${groupsize} \
                 --obsmaps \
                 --mem_count.enable \
                 --scan_map.disable \
                 --sim_atmosphere.cache_dir atm_cache_${telescope} \
                 --mapmaker.report_memory \
                 --save_hdf5.enable \
                 --save_hdf5.volume $outdir \
                 >> ${logfile} 2>&1 &
            let running++
            [[ $running -eq $nnode ]] && break
            #date >> ${logfile}
            #echo "Done with $logfile at" `date`
        else
            echo "$logfile exists"
        fi
    done
    [[ $running -eq $nnode ]] && break
done

echo "Waiting for $running jobs to complete at" `date`

wait

echo "Jobs completed at" `date`
