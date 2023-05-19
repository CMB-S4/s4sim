#!/bin/bash
#SBATCH --partition=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=64
#SBATCH --job-name=DC1_trial
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
export TOAST_LOGLEVEL=DEBUG

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=64
let groupsize=nnode_group*ntask_node

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "        groupsize = ${groupsize}"

logfile=trial.splat.log

if [[ ! -e $logfile ]]; then
    echo "Writing $logfile at" `date`
    date > ${logfile}
    srun -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
         --config trial.toml scanning_splat.toml atmosphere_pole.toml reduce_splat.toml \
         --focalplane focalplane_LAT2_SPLAT_f090.h5 \
         --schedule schedule.splat.short.txt \
         --out trial_outputs_splat \
         --job_group_size ${groupsize} \
         --obsmaps \
         --mem_count.enable \
         --mapmaker.report_memory \
         --scan_map.file cmb.splat.f090.h5 \
         --sim_atmosphere.cache_dir trial_atm_cache_splat \
         >> ${logfile} 2>&1
    date >> ${logfile}
    echo "Done with $logfile at" `date`
else
    echo "$logfile exists"
fi

#         --sim_noise.serial
