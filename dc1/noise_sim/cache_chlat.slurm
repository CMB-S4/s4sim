#!/bin/bash
#SBATCH --partition=regular
#SBATCH --time=04:00:00
#SBATCH --nodes=256
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
export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export TOAST_FUNCTIME=1
#export TOAST_LOGLEVEL=DEBUG

export PYTHONPATH=/global/common/software/cmb/keskital/cmbenv-intel/lib/python3.7/site-packages:/global/common/software/cmb/cori/cmbenv-intel_20211021/cmbenv_aux/lib/python3.7/site-packages

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS
let groupsize=4

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"

telescope=chlat
site=chile
TELESCOPE=LAT0_CHLAT
schedule=../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned.txt

logdir=logs/cache/${telescope}
mkdir -p $logdir

logdir=logs/cache
outdir=outputs/cache/${telescope}
mkdir -p $outdir
logfile=$logdir/cache_${telescope}.log
if [[ ! -e $logfile ]]; then
    echo "Writing $logfile at" `date`
    date > ${logfile}

    srun -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
         --config common.toml scanning_${telescope}.toml atmosphere_${site}.toml reduce_${telescope}.toml \
         --focalplane dummy_focalplane.440Hz.h5 \
         --schedule $schedule \
         --out $outdir \
         --job_group_size ${groupsize} \
         --mem_count.enable \
         --scan_map.disable \
         --sim_noise.disable \
         --sim_atmosphere.enable \
         --sim_atmosphere.cache_dir atm_cache_${telescope} \
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
