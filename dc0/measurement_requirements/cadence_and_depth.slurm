#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# OOM: 12
# Success: 
#SBATCH --cpus-per-task=24
#SBATCH --job-name=CMBS4_DC1_coadd
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
export OMP_NUM_THREADS=12

echo "Running with"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"

python cadence_and_depth.py
