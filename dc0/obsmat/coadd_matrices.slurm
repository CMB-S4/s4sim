#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --job-name=CMBS4_DC0_coadd_matrices
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
export OMP_NUM_THREADS=16
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
indirs=outputs/*/${band}/*
obsmats=""
for indir in $indirs; do
    obs=$(basename $indir)
    obsmat=$indir/filterbin_${obs}_noiseweighted_obs_matrix.npz
    # Get the multiplicity from the schedule file
    let n=`tail -n 1 split_schedule/${obs}.txt | awk '{print $10}'`
    if [[ -e $obsmat ]]; then
        echo "$(date) : Found: $obsmat $n"
        obsmats+=" ${obsmat}+${n}"
    else
        echo "$(date) : Not found: $obsmat"
        exit
    fi
done

outdir="obsmats/full/${band}"
mkdir -p $outdir
obsmat_full="$outdir/obsmat_${band}.npz"
invcov_full="$outdir/invcov_${band}.fits"
cov_full="$outdir/cov_${band}.fits"

logfile=${obsmat_full/.npz/.log}
echo "$(date) : Writing ${logfile} ..."
toast_obsmatrix_coadd \
    --outmatrix $obsmat_full \
    --invcov $invcov_full \
    --cov $cov_full \
    --rcond_limit 1e-3 \
    --double_precision \
    $obsmats \
    >& ${logfile}
echo "$(date) : ...Done!"
