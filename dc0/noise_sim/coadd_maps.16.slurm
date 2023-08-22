#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# OOM: 24
# Success: 48, 36, 30, 26
#SBATCH --cpus-per-task=26
#SBATCH --job-name=CMBS4_DC0_coadd_noise_16
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
export OMP_NUM_THREADS=13
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo "Running with"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"

indir=/global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk
let nsplit=32

# Random wait time to reduce clashes
# sleep $((RANDOM % 15))

#for TELESCOPE in LAT0_CHLAT; do
for TELESCOPE in SAT3_SAT; do
    case $TELESCOPE in
        LAT0_CHLAT)
            telescope=chlat
            bands=(f030 f040 f090 f150 f220 f280)
            ;;
        LAT2_SPLAT)
            telescope=splat
            bands=(f030 f040 f090 f150 f220 f280)
            ;;
        SAT1_SAT)
            telescope=spsat
            bands=(f095 f155 f220 f280)
            ;;
        SAT2_SAT)
            telescope=spsat
            bands=(f085 f095 f145 f155 f220 f280)
            ;;
        SAT3_SAT)
            telescope=spsat
            bands=(f030 f040 f085 f145)
            ;;
        *)
            echo "Unknown TELESCOPE: $TELESCOPE"
            ;;
    esac

    outdir=${indir}/coadd/${TELESCOPE}
    for band in ${bands[*]}; do
        # Check that all the input maps are in place
        for isplit in `seq 1 $nsplit`; do
            echo ${TELESCOPE} ${band}GHz $isplit
            splitroot=`printf "%03iof%03i" ${isplit} ${nsplit}`
            outroot=$outdir/coadd_${TELESCOPE}_${band}_${splitroot}
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
                outroot=$outdir/coadd_${TELESCOPE}_${band}_${splitroot}
                outmap=${outroot}_map.fits
                if [[ -e $outmap ]]; then
                    echo "$outmap already exists, skipping..."
                    continue
                fi
                logdir=coadd_logs
                mkdir -p $logdir
                logfile=${logdir}/coadd_${TELESCOPE}_${band}_${splitroot}.log
                if [[ -e $logfile ]]; then
                    echo "$logfile already exists, skipping..."
                    continue
                fi
                let isplit_in1=$isplit
                let isplit_in2=$isplit+$nsplit_out
                splitroot1=`printf "%03iof%03i" ${isplit_in1} ${nsplit_in}`
                splitroot2=`printf "%03iof%03i" ${isplit_in2} ${nsplit_in}`
                inroot1=$outdir/coadd_${TELESCOPE}_${band}_${splitroot1}
                inroot2=$outdir/coadd_${TELESCOPE}_${band}_${splitroot2}
                echo "Writing $logfile"
                date > $logfile
                toast_healpix_coadd \
                     --outmap ${outmap} \
                     --rcond ${outroot}_rcond.fits \
                     --rcond_limit 1e-3 \
                     --invcov ${outroot}_invcov.fits \
                     --cov ${outroot}_cov.fits \
                     ${inroot1}_map.fits \
                     ${inroot2}_map.fits \
                     >> $logfile 2>&1
            done
            let nsplit_in=$nsplit_out
        done
    done
done
