#!/bin/bash
#SBATCH --partition=regular
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --job-name=DC1_cache
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
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread
export TOAST_FUNCTIME=1
export TOAST_LOGLEVEL=DEBUG
#export KMP_AFFINITY=none,verbose

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS
let groupsize=1

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"

telescope=chlat
site=chile
TELESCOPE=LAT0_CHLAT
band=f220

fnames=split_schedules_32/${telescope}/split_schedule_????.txt

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

echo "Looking for schedule at" `date`

logdir=logs/cache/${telescope}
mkdir -p $logdir

# srun -n $ntask -c $ncore --cpu_bind=cores check-hybrid.intel.cori > thread_affinity.log

#export OMP_NUM_THREADS=1

for schedule in $fnames; do
    rootname=`basename $schedule .txt`
    rootname=${telescope}_${rootname}
    logfile=$logdir/cache_${rootname}.log
    outdir=outputs/cache/${rootname}
    if [[ ! -e $logfile ]]; then
        echo "Writing $logfile at" `date`
        date > ${logfile}

        #srun -n $ntask -c $ncore --cpu_bind=cores ../openmp_test/openmp_test >& openmp_test.log

        #srun -n $ntask -c $ncore --cpu_bind=cores python3 ../openmp_test/openmp_test.py >& openmp_test.log

        #exit
        
        srun -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
             --config common.toml scanning_${telescope}.toml atmosphere_${site}.toml reduce_${telescope}.toml \
             --focalplane ../focalplanes/focalplane_${TELESCOPE}_${band}.h5 --thinfp 1000000 --sample_rate 11 \
             --telescope $TELESCOPE \
             --schedule $schedule \
             --out $outdir \
             --job_group_size ${groupsize} \
             --mem_count.enable \
             --scan_map.disable \
             --sim_noise.disable \
             --sim_atmosphere.enable \
             --sim_atmosphere.cache_dir atm_cache_${TELESCOPE} \
             --sim_atmosphere.cache_only \
             --raw_statistics.disable \
             --filtered_statistics.disable \
             --mapmaker.report_memory \
             --save_hdf5.disable \
             >> ${logfile} 2>&1
        date >> ${logfile}
        echo "Done with $logfile at" `date`
        exit
    else
        echo "$logfile exists"
    fi
done

#             --sim_atmosphere.overwrite_cache
