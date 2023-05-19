#!/bin/bash
#SBATCH --partition=debug
#SBATCH --time=00:30:00
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
let nnode_group=4
let groupsize=nnode_group*ntask_node

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "        groupsize = ${groupsize}"

logfile=trial.spsat.log

if [[ ! -e $logfile ]]; then
    echo "Writing $logfile at" `date`
    date > ${logfile}
    srun -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
         --config trial.toml scanning_spsat.toml atmosphere_pole.toml reduce_spsat.toml \
         --focalplane focalplane_SAT0_SAT_f085.h5 \
         --schedule schedule.spsat.txt \
         --out trial_outputs_spsat \
         --job_group_size ${groupsize} \
         --obsmaps \
         --mem_count.enable \
         --mapmaker.report_memory \
         --scan_map.file cmb.spsat.f085.h5 \
         --pixels_radec.nside 512 \
         --sim_atmosphere.field_of_view "30 deg" \
         --sim_atmosphere.cache_dir trial_atm_cache_spsat \
         >> ${logfile} 2>&1
    date >> ${logfile}
    echo "Done with $logfile at" `date`
else
    echo "$logfile exists"
fi

#         --sim_noise.serial
#         --sim_atmosphere.overwrite_cache \
