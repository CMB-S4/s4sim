#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=08:00:00
#SBATCH --nodes=4
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

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim
indir1=${indir_root}/outputs_rk
#indir1=outputs
indir2=${indir_root}/outputs_DUMMY1
indir3=${indir_root}/outputs_DUMMY2
let nsplit=32

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    for band in f030 f040 f090 f150 f220 f280; do
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
        for isplit in `seq 1 $nsplit`; do
            echo ${telescope} ${band}GHz $isplit
            splitroot=`printf "_%03iof%03i" ${isplit} ${nsplit}`
            outroot=$outdir/coadd_${telescope}_${band}${splitroot}
            outmap=${outroot}_map.fits
            if [[ -e $outmap ]]; then
                echo "$outmap already exists, skipping..."
                continue
            fi
            logfile=coadd_${telescope}_${band}${splitroot}.log
            date > $logfile
            input_maps=""
            let ntotal=0
            let nfail=0
            let nfound=0
            fname_maps="coadd_maps_${band}${splitroot}.txt"
            rm -f $fname_maps
            for schedule in $schedule_dir/chlat/*txt; do
                if [[ $nsplit -gt 1 ]]; then
                    # Use the Julian date for splitting
                    let testsplit=`awk "{if (NR == 4) print (int(\\$5) % $nsplit) + 1}" $schedule`
                    [[ ! $testsplit -eq $isplit ]] && continue
                fi

                obs=`basename --suffix=.txt $schedule`
                fname1="${indir1}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
                fname2="${indir2}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
                fname3="${indir3}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
                for fname in $fname1 $fname2 $fname3 FAILED; do
                    [[ -e $fname ]] && break
                done
                if [[ $fname == FAILED ]]; then
                    echo "No input map for $band $obs" >> $logfile
                    let nfail++
                else
                    echo $fname >> $fname_maps
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

        # Stop at this level for now.  There has to be a more efficient way to
        # combine just two maps
        continue
        
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
                date > $logfile
                echo "Writing $logfile"
                date >> $logfile
                srun -N 1 -n 2 --cpu_bind=cores toast_healpix_coadd \
                     --outmap ${outmap} \
                     --rcond ${outroot}_rcond.fits \
                     --rcond_limit 1e-3 \
                     --invcov ${outroot}_invcov.fits \
                     --cov ${outroot}_cov.fits \
                     ${inroot1}_map.fits \
                     ${inroot2}_map.fits \
                     >> $logfile 2>&1 &
            done
            wait
            let nsplit_in=$nsplit_out
        done
    done
done
