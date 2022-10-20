#!/bin/bash
#SBATCH --partition=regular
#SBATCH --time=12:00:00
#SBATCH --nodes=128
#SBATCH --job-name=find_duplicates
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
export HDF5_USE_FILE_LOCKING=FALSE

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"

logfile=find_duplicates_in_splits.log
date >> $logfile
srun -n $ntask -c $ncore --cpu_bind=cores \
     python3 find_duplicates_in_splits.py \
     >> $logfile 2>&1
date >> $logfile
