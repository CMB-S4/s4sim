#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# OOM: 24
# Success: 48, 36, 30, 26
#SBATCH --cpus-per-task=26
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
export OMP_NUM_THREADS=13
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo "Running with"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim
indir1=${indir_root}/outputs_rk
#indir1=outputs
indir2=${indir_root}/outputs_DUMMY1
indir3=${indir_root}/outputs_DUMMY2
let nsplit=32

# Random wait time to reduce clashes
# sleep $((RANDOM % 15))

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    for band in f030 f040 f090 f150 f220 f280; do
    #for band in f280; do
        case $band in
            f030|f040|f090|f150)
                schedule_dir=split_schedules_1_upto3mm_with_break
                ;;
            f220|f280)
                schedule_dir=split_schedules_1_upto2mm_with_break
                ;;
            *)
                echo "Unknown frequency band: $band"
                exit
                ;;
        esac

        # Check that all the input maps are in place
        for isplit in `seq 1 $nsplit`; do
            echo ${telescope} ${band}GHz $isplit
            splitroot=`printf "_%03iof%03i" ${isplit} ${nsplit}`
            outroot=$outdir/coadd_${telescope}_${band}${splitroot}
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
                splitroot=`printf "_%03iof%03i" ${isplit} ${nsplit_out}`
                outroot=$outdir/coadd_${telescope}_${band}${splitroot}
                outmap=${outroot}_map.fits
                if [[ -e $outmap ]]; then
                    echo "$outmap already exists, skipping..."
                    continue
                fi
                let isplit_in1=$isplit
                let isplit_in2=$isplit+$nsplit_out
                splitroot1=`printf "_%03iof%03i" ${isplit_in1} ${nsplit_in}`
                splitroot2=`printf "_%03iof%03i" ${isplit_in2} ${nsplit_in}`
                inroot1=$outdir/coadd_${telescope}_${band}${splitroot1}
                inroot2=$outdir/coadd_${telescope}_${band}${splitroot2}
                logfile=coadd_${telescope}_${band}${splitroot}.log
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
