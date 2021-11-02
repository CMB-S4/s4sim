#!/bin/bash
#SBATCH --partition=regular
#SBATCH --time=04:00:00
#SBATCH --nodes=32
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

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=32
let groupsize=nnode_group*ntask_node

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "        groupsize = ${groupsize}"

logfile=trial.chlat.log

if [[ ! -e $logfile ]]; then
    echo "Writing $logfile at" `date`
    date > ${logfile}
    srun -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
         --config common.toml scanning_chlat.toml atmosphere_chile.toml reduce_chlat.toml \
         --focalplane focalplane_LAT0_CHLAT_f150.h5 \
         --schedule schedule.short.txt \
         --out trial_outputs_chlat \
         --job_group_size ${groupsize} \
         --obsmaps \
         --mem_count.enable \
         --mapmaker.report_memory \
         >> ${logfile} 2>&1
    date >> ${logfile}
    echo "Done with $logfile at" `date`
else
    echo "$logfile exists"
fi

#         --sim_noise.serial
