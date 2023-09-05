#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# OOM:
# Success:
#SBATCH --cpus-per-task=8
#SBATCH --job-name=CMBS4_DC0_wrap_up_obs
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
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo "Running with"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"

#python3 wrap_up_obs_maps.py 025
#python3 wrap_up_obs_maps.py 040
#python3 wrap_up_obs_maps.py 090
#python3 wrap_up_obs_maps.py 150
#python3 wrap_up_obs_maps.py 230
python3 wrap_up_obs_maps.py 280
