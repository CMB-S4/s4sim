#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# OOM: 24
# Success: 48, 36, 30, 26
#SBATCH --cpus-per-task=32
#SBATCH --job-name=CMBS4_DC0_coadd_16
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107

# Run this script *after* the 32-way splits have been produced with
# coadd_maps.slurm
# This script generates the 16/8/4/2/1 splits

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

echo "Running with"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"

indir=/global/cfs/cdirs/cmbs4/dc/dc0/staging/multimap_sim/outputs_rk
let nsplit=32

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

for telescope in LAT0_CHLAT; do
    outdir=$indir/coadd/${telescope}
    #for band in f030 f040 f090 f150 f220 f280; do
    for band in f280; do
        # Check that all the input maps are in place
        for flavor in unlensed_cmb cmb_lensing foreground; do
            for isplit in `seq 1 $nsplit`; do
                echo ${telescope} ${band}GHz $isplit
                splitroot=`printf "%03iof%03i" ${isplit} ${nsplit}`
                outroot=$outdir/coadd_${telescope}_${band}_${flavor}_${splitroot}
                outmap=${outroot}_map.fits
                if [[ ! -e $outmap ]]; then
                    echo "$outmap does not exist!"
                    exit
                fi
            done

            # Co-add maps
            let nsplit_in=$nsplit
            while [[ $nsplit_in -gt 1 ]]; do
                let nsplit_out=$nsplit_in/2
                for isplit in `seq 1 $nsplit_out`; do
                    splitroot=`printf "%03iof%03i" ${isplit} ${nsplit_out}`
                    outroot=$outdir/coadd_${telescope}_${band}_${flavor}_${splitroot}
                    outmap=${outroot}_map.fits
                    if [[ -e $outmap ]]; then
                        echo "$outmap already exists, skipping..."
                        continue
                    fi
                    logdir=coadd_logs
                    mkdir -p $logdir
                    logfile=${logdir}/coadd_${telescope}_${band}_${flavor}_${splitroot}.log
                    if [[ -e $logfile ]]; then
                        echo "$logfile already exists, skipping..."
                        continue
                    fi
                    let isplit_in1=$isplit
                    let isplit_in2=$isplit+$nsplit_out
                    splitroot1=`printf "%03iof%03i" ${isplit_in1} ${nsplit_in}`
                    splitroot2=`printf "%03iof%03i" ${isplit_in2} ${nsplit_in}`
                    inroot1=$outdir/coadd_${telescope}_${band}_${flavor}_${splitroot1}
                    inroot2=$outdir/coadd_${telescope}_${band}_${flavor}_${splitroot2}
                    # Link the inverse covariance from noise_sim
                    invcov1_in=${outdir/multimap_sim/noise_sim}/coadd_${telescope}_${band}_${splitroot1}_invcov.fits
                    invcov2_in=${outdir/multimap_sim/noise_sim}/coadd_${telescope}_${band}_${splitroot2}_invcov.fits
                    invcov1_out=${inroot1}_invcov.fits
                    invcov2_out=${inroot2}_invcov.fits
                    ln -s $invcov1_in $invcov1_out
                    ln -s $invcov2_in $invcov2_out
                    echo "Writing $logfile"
                    date > $logfile
                    toast_healpix_coadd \
                        --outmap ${outmap} \
                        --rcond_limit 1e-3 \
                        ${inroot1}_map.fits \
                        ${inroot2}_map.fits \
                        >> $logfile 2>&1
                done # isplit
                let nsplit_in=$nsplit_out
            done # nsplit
        done # flavor
    done # band
done # telescope
