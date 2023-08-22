#!/bin/bash
#SBATCH --qos=preempt
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --job-name=CMBS4_DC0_coadd_noise_32
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107

# This script produces the 32-way split maps from
# single observation maps

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
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"

indir=/global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk
let nsplit=32

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

#for TELESCOPE in LAT0_CHLAT; do
for TELESCOPE in SAT1_SAT SAT2_SAT SAT3_SAT; do
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
            # bands=(f030 f040)
            ;;
        *)
            echo "Unknown TELESCOPE: $TELESCOPE"
            ;;
    esac

    outdir=${indir}/coadd/${TELESCOPE}
    # Determine the appropriate scaling to match nominal sensitivity
    case $TELESCOPE in
        LAT0_CHLAT)
            # Scale the noise maps to match 7-year 2-LAT sensitivity
            scale=$(echo "1 / sqrt(2 * 7)" | bc -l)
            ;;
        *)
            # Scale the noise maps to match 7-year sensitivity
            scale=$(echo "1 / sqrt(1 * 7)" | bc -l)
            ;;
    esac
    for band in ${bands[*]}; do
        if [[ $TELESCOPE == "LAT0_CHLAT" ]]; then
            # For CHLAT we simulated observations across the year but
            # only want to co-add the good weather nominal season
            case $band in
                f030|f040|f090|f150)
                    schedule_dir=../split_schedules_1_upto3mm_with_break
                    ;;
                f220|f280)
                    schedule_dir=../split_schedules_1_upto2mm_with_break
                    ;;
                *)
                    echo "Unknown frequency band: $band"
                    exit
                    ;;
            esac
        else
            # For SPSAT and SPLAT we only simulated the nominal season and all
            # simulated observations make the cut
            schedule_dir=../split_schedules_1_upto2mm
        fi
        for isplit in `seq 1 $nsplit`; do
            echo ${TELESCOPE} ${band}GHz $isplit
            splitroot=`printf "%03iof%03i" ${isplit} ${nsplit}`
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
            date > $logfile
            input_maps=""
            let ntotal=0
            let nfail=0
            let nfound=0
            listdir=map_lists
            mkdir -p $listdir
            fname_maps=${listdir}/coadd_maps_${band}_${splitroot}.txt
            rm -f $fname_maps
            for schedule in $schedule_dir/$telescope/*txt; do
                if [[ $nsplit -gt 1 ]]; then
                    # Use the Julian date for splitting.
                    let testsplit=`awk "{if (NR == 4) print (int(\\$5) % $nsplit) + 1}" $schedule`
                    [[ ! $testsplit -eq $isplit ]] && continue
                fi

                obs=`basename --suffix=.txt $schedule`
                if [[ $TELESCOPE == SAT* ]]; then
                    fname="${indir}/${TELESCOPE}/${band}/${obs}/filterbin_${obs}_noiseweighted_filtered_map.h5"
                    # [[ ! -e $fname ]] && fname="${indir}/${TELESCOPE}/${band}/${obs}/filterbin_noiseweighted_filtered_map.h5"
                else
                    fname="${indir}/${TELESCOPE}/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
                fi
                if [[ ! -e $fname ]]; then
                    echo "Not found: $fname"
                    echo "Not found: $fname" >> $logfile
                    exit
                    let nfail++
                else
                    echo $fname >> $fname_maps
                    let nfound++
                fi
                let ntotal++
            done
            if [[ ! -e $fname_maps ]]; then
                echo "ERROR: $fname_maps was not created."
                continue
            fi
            echo "Found a total of ${nfound} / ${ntotal} maps. ${nfail} maps were missing." >> $logfile

            mkdir -p $outdir
            echo "Writing $logfile"
            date >> $logfile
            srun -n $ntask -c $ncore --cpu_bind=cores toast_healpix_coadd \
                 --scale ${scale} \
                 --outmap ${outmap} \
                 --rcond ${outroot}_rcond.fits \
                 --rcond_limit 1e-3 \
                 --invcov ${outroot}_invcov.fits \
                 --cov ${outroot}_cov.fits \
                 $fname_maps \
                 >> $logfile 2>&1
            date >> $logfile
        done
    done
done
