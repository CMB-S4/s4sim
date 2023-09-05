#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# OOM: 8, 16, 24, 32
# Success:
#SBATCH --cpus-per-task=64
#SBATCH --job-name=CMBS4_DC0_combine_matrices
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

echo "$(date) : Running with OMP_NUM_THREADS = ${OMP_NUM_THREADS}"

#band=f030
#band=f040
#band=f085
#band=f095
#band=f145
#band=f155
#band=f220
band=f280
indirs=outputs/*/$band/*
obsmats=""
for indir in $indirs; do
    obs=$(basename $indir)
    obsmat=$indir/filterbin_${obs}_noiseweighted_obs_matrix.npz
    obsmats+=" $obsmat"
    if [[ -e $obsmat ]]; then
        echo "$(date) : Found: $obsmat, skipping..."
        continue
    fi
    logfile=${obsmat/.npz/.log}
    if [[ -e $logfile ]]; then
        echo "$(date) : $logfile exists, skipping..."
        continue
    fi
    echo "$(date) : Writing ${logfile} ..."
    toast_obsmatrix_combine ${obsmat/.npz/} >& ${logfile}
    echo "$(date) : ... Done!"
done
echo "$(date) : All done."
