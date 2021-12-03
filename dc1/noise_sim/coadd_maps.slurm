#!/bin/bash
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --job-name=coadd_maps
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

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    for band_dir in ${indir_root}/${telescope}/f???; do
        band=`basename $band_dir`
        echo ${telescope} ${band}GHz
        indir=${indir_root}/${telescope}/${band}
        mkdir -p $outdir
        outroot=$outdir/coadd_${telescope}_${band}
        logfile=coadd_${telescope}_${band}.log
        echo "Writing $logfile"
        srun -n $ntask -c $ncore --cpu_bind=cores toast_healpix_coadd \
            --outmap ${outroot}_map.fits \
            --rcond ${outroot}_rcond.fits \
            --rcond_limit 1e-3 \
            --invcov ${outroot}_invcov.fits \
            --cov ${outroot}_cov.fits \
            ${indir}/${telescope}_split_schedule_00??/mapmaker_*_noiseweighted_map.h5 \
            >& $logfile
    done
done
