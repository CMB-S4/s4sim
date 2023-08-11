#!/bin/bash
#SBATCH --qos=preempt
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --job-name=CMBS4_DC0_coadd_32
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

indir=/global/cfs/cdirs/cmbs4/dc/dc0/staging/multimap_sim/outputs_rk
let nsplit=32

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

for telescope in LAT0_CHLAT; do
    case $telescope in
        LAT0_CHLAT)
            tele=chlat
            ;;
        LAT2_SPLAT)
            tele=splat
            ;;
        *)
            echo "Unknown telescope: $telescope"
            ;;
    esac
    outdir=$indir/coadd/${telescope}
    #for band in f030 f040 f090 f150 f220 f280; do
    for band in f030 f040; do
        if [[ $telescope == "LAT0_CHLAT" ]]; then
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
            echo ${telescope} ${band}GHz $isplit
            splitroot=`printf "%03iof%03i" ${isplit} ${nsplit}`
            for flavor in unlensed_cmb cmb_lensing foreground; do
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
                date > $logfile
                input_maps=""
                let ntotal=0
                let nfail=0
                let nfound=0
                listdir=map_lists
                mkdir -p $listdir
                fname_maps=${listdir}/coadd_maps_${band}_${flavor}_${splitroot}.txt
                rm -f $fname_maps
                for schedule in $schedule_dir/$tele/*txt; do
                    if [[ $nsplit -gt 1 ]]; then
                        # Use the Julian date for splitting
                        let testsplit=`awk "{if (NR == 4) print (int(\\$5) % $nsplit) + 1}" $schedule`
                        [[ ! $testsplit -eq $isplit ]] && continue
                    fi

                    obs=`basename --suffix=.txt $schedule`
                    fname="${indir}/${telescope}/${band}/${obs}/mapmaker_${obs}_${flavor}_noiseweighted_map.h5"
                    if [[ ! -e $fname ]]; then
                        echo "Not found: $fname" >> $logfile
                        let nfail++
                    else
                        echo $fname >> $fname_maps
                        # Ensure we have a symbolic link in place to the
                        # inverse white noise covariance
                        fname_invcov_out=${fname/noiseweighted_map/invcov}
                        if [[ ! -f $fname_invcov_out ]]; then
                            rm -f $fname_invcov_out
                            fname_invcov_in=${fname_invcov_out/_$flavor/}
                            fname_invcov_in=${fname_invcov_in/multimap_sim/noise_sim}
                            if [[ ! -e $fname_invcov_in ]]; then
                                echo "ERROR: missing invcov: $fname_invcov_in"
                                exit
                            fi
                            ln -s $fname_invcov_in $fname_invcov_out
                        fi
                        let nfound++
                    fi
                    let ntotal++
                done
                echo "Found a total of ${nfound} / ${ntotal} maps. ${nfail} maps were missing." >> $logfile

                mkdir -p $outdir
                echo "Writing $logfile"
                date >> $logfile
                srun -n $ntask -c $ncore --cpu_bind=cores toast_healpix_coadd \
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
done
